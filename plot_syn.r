suppressPackageStartupMessages(require(argparse))

ps <- ArgumentParser(description = 'Syntenic visualization using provided msa file')
ps$add_argument("fi", help="msa file")
ps$add_argument("fo", help="output prefix")

ps$add_argument("--nwk", help="nwk tree")
ps$add_argument("--gffdir", help="gff file dir")
ps$add_argument("--up", type='integer', help="upstream", default=5000)
ps$add_argument("--down", type='integer',  help="downstream", default=3000)
args <- ps$parse_args()

#setwd("/data/jinlab/Public/LYF_haplotype/NAM")
require(Biostrings)
library("dplyr")
require(tidyverse)
require("stringr")
library("purrr")
library("ggpubr")
library("ggplot2")
library("ggtree")
library("magrittr")
library("limma")
require(glue)


msa_tree <- function(nwk, seqs, aligner) {

  tree = read.tree(nwk)
  labs = with(subset(fortify(tree), isTip), label[order(y, decreasing=T)])
  nl = length(labs)
  ty = tibble(gt = labs, y = nl:1) %>%
    #separate(gt, c('i','gt'), extra='merge', sep='_') %>%
    mutate(lab = str_replace(gt, '^Zmays_', '')) %>%
    mutate(t = strsplit2(gt,"-")[,2]) %>%
    mutate(g = strsplit2(gt,"_")[,1]) %>%
    dplyr::select(gt, lab, y, t, g)
  #

  list(tree=tree, ty=ty)
  
  #ty2 = ty %>% select(taxa=lab, gt)
  #ggtree(tree) %<+% ty2 +
    #geom_tiplab(aes(label=gt), hjust=0, align=T) +
    #scale_x_continuous(expand=expansion(mult=c(.05,.005))) +
    #scale_y_continuous(expand=expansion(mult=c(.01,.05)))
    
  #}}}
}


make_syn  <- function(name1, name2, seqs, ht) {
  #{{{ extract synteny blocks and mismatch positions from alignment
  #x = unmasked(msa)
  #s1 = RemoveGaps(x[name1])[[1]]
  #s2 = RemoveGaps(x[name2])[[1]]
  name1 <- names(seqs)[grep(name1,names(seqs))]
  name2 <- names(seqs)[grep(name2,names(seqs))]
  n1st <- strsplit(strsplit(name1,split=":")[[1]][4], split="-")[[1]][1]
  n1ed <- strsplit(strsplit(name1,split=":")[[1]][4], split="-")[[1]][2]
  n2st <- strsplit(strsplit(name2,split=":")[[1]][4], split="-")[[1]][1]
  n2ed <- strsplit(strsplit(name2,split=":")[[1]][4], split="-")[[1]][2]
  
  
  
  stopifnot(name1 %in% names(seqs) & name2 %in% names(seqs))
  pw = pairwiseAlignment(seqs[[name1]], seqs[[name2]], type='global')
  toff = start(pattern(pw)); qoff = start(subject(pw))
  print(toff); print(qoff)
  y1 = as.character(pattern(pw))
  y2 = as.character(subject(pw))
  h0 = as.character(compareStrings(y2,y1))
  h = str_replace_all(h0, "[ATCGN?]", '*')
  h2=unlist(strsplit(h, split = ""))
  j = tibble(v = rle(h2)$values, aSize = rle(h2)$lengths) %>%
    mutate(aEnd=cumsum(aSize), aBeg = aEnd-aSize) %>%
    mutate(tSize = ifelse(v=='+', 0, aSize)) %>%
    mutate(qSize = ifelse(v=='-', 0, aSize)) %>%
    mutate(tEnd=cumsum(tSize), tBeg=tEnd-tSize) %>%
    mutate(qEnd=cumsum(qSize), qBeg=qEnd-qSize) %>%
    dplyr::select(v, aBeg,aEnd,aSize,tBeg,tEnd,tSize,qBeg,qEnd,qSize) %>%
    mutate(tBeg = toff+tBeg-1, tEnd=tBeg+tSize) %>%
    mutate(qBeg = qoff+qBeg-1, qEnd=qBeg+qSize)
  get_tpos <- function(v, aSize, tBeg, tEnd)
    ifelse(v=='+', list(rep(tBeg, aSize)), list((tBeg+1):tEnd))
  get_qpos <- function(v, aSize, qBeg, qEnd)
    ifelse(v=='-', list(rep(qBeg, aSize)), list((qBeg+1):qEnd))
  j2 = j %>%
    mutate(tpos = pmap(list(v, aSize, tBeg, tEnd), get_tpos)) %>%
    mutate(qpos = pmap(list(v, aSize, qBeg, qEnd), get_qpos))
  tposs = unlist(j2$tpos)
  qposs = unlist(j2$qpos)
  mm = str_locate_all(h0, "\\?")[[1]] %>% as_tibble() %>%
    mutate(qPos = qposs[start], tPos = tposs[start]) %>%
    dplyr::select(aPos=start,tPos,qPos)
  list(aln = j, mm = mm)
  j <- as.data.frame(j)
  #j$tBeg <- j$tBeg + as.numeric(n1st)
  #j$tEnd <- j$tEnd + as.numeric(n1st)
  #j$qBeg <- j$qBeg + as.numeric(n2st)
  #j$qEnd <- j$qEnd + as.numeric(n2st)
  t <- j[j$v=="*",]
  n = 1
  ty = ht-0.1
  qy = ht-0.9
  for (i in 1:nrow(t)) {
    #print (c(t[i,]$tBeg,t[i,]$tEnd,t[i,]$qEnd,t[i,]$qBeg))
    tb <- data.frame(pos=c(t[i,]$tBeg,t[i,]$tEnd,t[i,]$qEnd,t[i,]$qBeg), y=c(ty,ty,qy,qy), group=c(n,n,n,n))
    if (n==1){out <- tb}else{out <- rbind(out, tb)}
    n=n+1
  }
  
  p_syn = ggplot(out) +
    geom_polygon(aes(x=pos,y=y,group=group), fill='royalblue', alpha=.2,
                 size=0,color=NA) 
  out
  #}}}
}


make_syn2  <- function(name1, name2, seqs, ht) {
  #{{{ extract synteny blocks and mismatch positions from alignment
  #x = unmasked(msa)
  #s1 = RemoveGaps(x[name1])[[1]]
  #s2 = RemoveGaps(x[name2])[[1]]
  name1 <- names(seqs)[grep(name1,names(seqs))]
  name2 <- names(seqs)[grep(name2,names(seqs))]
  n1st <- strsplit(strsplit(name1,split=":")[[1]][4], split="-")[[1]][1]
  n1ed <- strsplit(strsplit(name1,split=":")[[1]][4], split="-")[[1]][2]
  n2st <- strsplit(strsplit(name2,split=":")[[1]][4], split="-")[[1]][1]
  n2ed <- strsplit(strsplit(name2,split=":")[[1]][4], split="-")[[1]][2]
  
  
  
  stopifnot(name1 %in% names(seqs) & name2 %in% names(seqs))
  #pw = pairwiseAlignment(seqs[[name1]], seqs[[name2]], type='global')
  #seqs = readDNAStringSet("/data/jinlab/Public/LYF_haplotype/NAM/lg1_4k_4k.2.fa")
  y1 <- as.character(seqs[[name1]])
  y2 <- as.character(seqs[[name2]])
  y1 <- str_replace_all(y1, "[ATCGN?]", '+')
  y2 <- str_replace_all(y2, "[ATCGN?]", '+')
  y1 <- unlist(strsplit(y1, split = ""))
  y2 <- unlist(strsplit(y2, split = ""))

  a1 = 0
  a2 = 0
  d = F
  o <- data.frame()
  for (i in 1:length(y1)) {
    if ((y1[i] == y2[i]) & (y1[i] == "-"))
    {#if(d==TRUE){ed1<-a1;ed2<-a2;print(paste(st1,st2,ed1,ed2,sep=","));d = F}else{};
      next;}
    else
    {if ((y1[i] == "+") & (y2[i] == "+"))
    {a1=a1+1; a2=a2+1; if(d==TRUE){}else{d=TRUE;st1<-a1;st2<-a2}}
      else
      {if(d==TRUE){ed1<-a1;ed2<-a2;o <- rbind(o, data.frame(st1,st2,ed1,ed2));d = F}else{};
        if (y1[i] == "+"){a1=a1+1}else{a2=a2+1}}
    }
  }
  if(d==TRUE){ed1<-a1;ed2<-a2;o <- rbind(o, data.frame(st1,st2,ed1,ed2));d = F}else{}
  
  n = 1
  ty = ht-0.2
  qy = ht-0.8
  for (i in 1:nrow(o)) {
    #print (c(t[i,]$tBeg,t[i,]$tEnd,t[i,]$qEnd,t[i,]$qBeg))
    tb <- data.frame(pos=c(o[i,]$st1,o[i,]$ed1,o[i,]$ed2,o[i,]$st2), y=c(ty,ty,qy,qy), group=c(n,n,n,n))
    if (n==1){out <- tb}else{out <- rbind(out, tb)}
    n=n+1
  }
  
  p_syn = ggplot(out) +
    geom_polygon(aes(x=pos,y=y,group=group), fill='royalblue', alpha=.2,
                 size=0,color=NA) 
  out
  #}}}
}

upstream<-args$up#5000
downstream<-args$down#3000
seqs = readDNAStringSet(args$fi)
msa <- msa_tree(args$nwk,seqs,"muscle")
outdir<-args$gffdir



n=1
for (e in 2:nrow(msa$ty)) {
  name1 <- msa$ty[e-1,]$lab
  name2 <- msa$ty[e,]$lab
  h <- as.numeric(msa$ty[e-1,]$y)
  o <- make_syn2(name1, name2, seqs, h)
  if (n==1)
    {out <- o}
  else{
    maxg <- max(out$group)
    o$group <- o$group + maxg
    out <- rbind(out, o)}
  n=0
}

seqlen <- c()
for (s in 1:length(seqs)) {
  seqlen <- c(seqlen, as.numeric(length(strsplit(as.character(seqs[[s]]),"[ATCGN?]")[[1]])))
}

xbrks=c(0,upstream,min(seqlen)); xlabs=c(paste0("-",upstream/1000,"k"),'TSS',paste0('+',downstream/1000,"k"))#; xmax = 4000
p_syn = ggplot(out) +
  geom_polygon(aes(x=pos,y=y,group=group), fill='moccasin', alpha=1,
               size=0,color=NA) +
  scale_x_continuous(breaks=xbrks, labels=xlabs,
                     expand=expansion(mult=c(.01,.01)), position='top') +
  scale_y_continuous(breaks=msa$ty$y, labels=msa$ty$t, expand=expansion(mult=c(.01,.01))) + theme_bw()+
  theme(legend.background = element_blank(),strip.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
        panel.border = element_blank(),plot.margin = unit(c(.2,.5,.2,0), "lines"))



col.intron='grey'; col.exon='orange'; col.cds='#0476D9'
arw = arrow(length=unit(.1,'cm'), angle=30, ends='last',type="open")


n=1
for (e in 1:nrow(msa$ty)) {
  et <- msa$ty[e,]$g
  gff <- read.table(paste0(outdir, et, ".gff"), sep = "\t", header = F)
  h <- as.numeric(msa$ty[e,]$y)
  gff <- gff[,c(1,3,4,5,7)]
  gff$h <- h
  if(gff$V7=="+"){
    ma <- min(gff[gff$V3=="gene",]$V4,gff[gff$V3=="gene",]$V5)
    gff$V4 <- gff$V4 - ma + upstream
    gff$V5 <- gff$V5 - ma + upstream
  }else{
    ma <- min(gff[gff$V3=="gene",]$V4,gff[gff$V3=="gene",]$V5)
    gff$V4 <- gff$V4 - ma + upstream
    gff$V5 <- gff$V5 - ma + upstream
  }
    
  #gff$type <- NA
  if(n==1){gffcom <- gff}else{gffcom <- rbind(gffcom, gff)}
  n=0
  
}


ht.exon <- 0.06; ht.cds <- 0.08
arw.beg = upstream
arw.end = arw.beg + max(seqlen)/40
p_aln = p_syn +
  geom_segment(data=gffcom[gffcom$V3=="gene",],aes(x=V4,xend=V5,y=h,yend=h),col=col.intron,size=.2) +
  geom_rect(data=gffcom[gffcom$V3=="exon",],aes(xmin=V4,xmax=V5,ymin=h-ht.exon,ymax=h+ht.exon),fill=col.exon,color=NA,alpha=1) +
  geom_rect(data=gffcom[gffcom$V3=="CDS",],aes(xmin=V4,xmax=V5,ymin=h-ht.cds,ymax=h+ht.cds),fill=col.cds,color=NA,alpha=1) +
  geom_segment(data=gffcom[gffcom$V3=="gene",],aes(x=arw.beg,xend=arw.end,y=h+ht.cds*1.1,yend=h+ht.cds*1.1),
               color='black', size=.2, arrow=arw) +
  geom_vline(xintercept = c(upstream), color="yellow")
  #geom_rect(xmin=upstream-upstream/3000,xmax=upstream+upstream/3000,ymin=-Inf,ymax=Inf, fill='yellow', alpha=.2)


p_tree = ggtree(msa$tree)+
  #geom_tiplab(aes(label=gt), hjust=0, align=T) +
  scale_x_continuous(expand=expansion(mult=c(.05,.005))) +
  scale_y_continuous(expand=expansion(mult=c(.01,.05))) + #theme_bw()+
  theme(plot.margin = unit(c(.2,.2,.2,0), "lines"))

p = ggarrange(p_tree, p_aln, nrow=1, ncol=2, widths=c(1,5))
p

p %>% ggexport(filename=glue(paste0(args$fo, ".pdf")), width=7, height=6)

#name1 <- c("Ki3")
#name2 <- c("CML228")
#o <- make_syn(name1, name2, seqs)

