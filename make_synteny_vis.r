#options(repr.plot.width = 6, repr.plot.height = 10)
# 导入argparse包
library(argparse)
library("ggplot2")
# 创建参数解析器
parser <- ArgumentParser(description = "Pangenome Synteny visualization")

# 添加命令行参数

parser$add_argument("--info", help = "Genomic information file", required = TRUE)
parser$add_argument("--upstream", type = "numeric", help = "Upstream (default=10000)", default=10000)
parser$add_argument("--downstream", type = "numeric", help = "Downstream (default=10000)", default=10000)
parser$add_argument("--zoomin", type = "numeric", help = "Zoomin (default=0)", default=0)
parser$add_argument("--out", help = "output file prefix", required = TRUE)
parser$add_argument("--width", type = "numeric", help = "Figure width (default=6)", default=6)
parser$add_argument("--height", type = "numeric", help = "Figure height (default=3)", default=3)
#parser$add_argument("--bedtools", type = "numeric", help = "Figure height", default=10)
#parser$add_argument("--height", type = "numeric", help = "Figure height", default=10)
parser$add_argument("--snp", action = "store_true", help = "Flag to indicate if the snp/indels should be showed", default=TRUE)
parser$add_argument("--run", action = "store_true", help = "Flag to indicate if the mummer should run", default=TRUE)
#parser$add_argument("--plotout", action = "store_true", help = "Flag to indicate if the script should plot")

# 解析命令行参数
args <- parser$parse_args()

# 使用解析后的参数
info <- args$info
zoomin <- args$zoomin
upstream <- args$upstream
downstream <- args$downstream
snpindel <- args$snp
run <- args$run
prefix <- args$out
hh <- args$height
ww <- args$width

# 打印参数信息
cat("info:", info, "\n")
cat("upstream:", upstream, "\n")
cat("downstream:", downstream, "\n")
cat("zoomin:", zoomin, "\n")
cat("run:", run, "\n")
cat("output:", prefix, "\n")

info <- read.table(info, header = T)

n = 1
block_len = 0.40
exon_h = 0.06

nn = 1.5

info

for (i in seq_len(nrow(info))){
    line <- info[i,]$Region
    #print (i)
    prefix <- info[i,]$Prefix
    gff <- info[i,]$GFF
    
    strand = info[i,]$Strand
    t <- strsplit(line,":")
    chro = t[[1]][1]
    start = as.numeric(strsplit(t[[1]][2],"-")[[1]][1]) - upstream
    end = as.numeric(strsplit(t[[1]][2],"-")[[1]][2]) + downstream
    fa <- paste0(prefix, "_", chro, "_", start, "_", end, ".fasta")
    cmd1 <- paste0("echo -e '", chro, "\t", start, "\t", end, "' | bedtools getfasta -bed - -fi ",
                 info[i,]$Genome, " -fo ", fa)
    if (run == TRUE){system(cmd1)}else{} ###########################################################################################################################################################
    
    if (strand == FALSE){
        cmdd <- paste0("~/software/seqkit seq -t dna -rp ", fa, " > tmp.fa && mv tmp.fa ", fa)
        print (cmdd)
        if (run == TRUE){system(cmdd)}else{}
        #break
    }else{}
    #print (fa)
    
    cmd2 <- paste0("echo -e '", chro, "\t", start, "\t", end, "' | bedtools intersect -a - -b ",
                 gff, " -wo -nonamecheck > ", paste0(prefix, ".gff3"))
    print (cmd2)
    if (run == TRUE){system(cmd2)}else{} ###########################################################################################################################################################
    
    
    gf <- data.frame()
    if (file.info(paste0(prefix, ".gff3"))$size != 0){
        gf <- read.table(paste0(prefix, ".gff3"), header = F, sep = "\t")
        gf$prefix <- prefix
        gf <- gf[(gf$V6=="gene")|(gf$V6=="exon"),]
        gf$exon_h1 <- nn - exon_h; gf$exon_h2 <- nn + exon_h; gf$gene_h <- nn
        if (strand == FALSE){
        #gf[, 10] <- ifelse(gf[, 10] == "-", "+", "-")
        gf$V8 <- gf$V3 - gf$V8 + gf$V2
        gf$V7 <- gf$V3 - gf$V7 + gf$V2}
        #print (head(gf))
    }else{
        TRUE
        print (prefix)
    }
    #print (gf)
    
    nn = nn -1
    
    if (n == 1){
        l_fa <- fa; l_chro <- chro; l_start <- start; l_end <- end; l_pre <- prefix
        n = 0
        
        gfftotal <- gf
    } else if (n == 0){
        out <- data.frame(prefix1 = l_pre, prefix2 = prefix, fa1 = l_fa, fa2 = fa, chr1 = l_chro, chr2 = chro, start1 = l_start, start2 = start, end1 = l_end, end2 = end)
        l_fa <- fa; l_chro <- chro; l_start <- start; l_end <- end; l_pre <- prefix
        n = -1
        gfftotal <- rbind(gfftotal, gf)
    } else {
        out <- rbind(out, data.frame(prefix1 = l_pre, prefix2 = prefix, fa1 = l_fa, fa2 = fa, chr1 = l_chro, chr2 = chro, start1 = l_start, start2 = start, end1 = l_end, end2 = end))
        l_fa <- fa; l_chro <- chro; l_start <- start; l_end <- end; l_pre <- prefix
        gfftotal <- rbind(gfftotal, gf)
    }
    #print (end)
    #system2()
}

gfftotal$arrow <- ""
gfftotal[gfftotal$V10=="+",]$arrow <- "last"
gfftotal[gfftotal$V10=="-",]$arrow <- "first"
#out
out$len1 <- as.numeric(out$end1) - as.numeric(out$start1)
out$len2 <- as.numeric(out$end2) - as.numeric(out$start2)

m <- max(c(max(out$len1), c(max(out$len2))))
#gfftotal
#break

#block_len = 0.45

n = 1
nn = 1
for (i in seq_len(nrow(out))){
    l <- out[i,]
    prefix = paste0(l$prefix1, "_", l$prefix2)
    cmd <- paste0("dnadiff -p ",  prefix, " ", l$fa1, " ", l$fa2)
    
    #print (cmd)
    if (run == TRUE){system(cmd)}else{print("1")} ###################################################################################################################################################
    
    
    mcoords <- read.table(paste0(prefix, ".mcoords"), header = F)
    
    if (file.info(paste0(prefix, ".snps"))$size == 0){
        snps <- data.frame()
    }else{
        snps <- read.table(paste0(prefix, ".snps"), header = F)
    }
    
    
    
    if (n == 1){
        mcoords$prefix <- paste0(l$prefix1, l$prefix2)
        mcoords$n <- n
        mt <- mcoords
    }else{
        mcoords$prefix <- paste0(l$prefix1, l$prefix2)
        mcoords$n <- n
        mt <- rbind(mt, mcoords)
    }
    
    n = n - 1
    
    
    
    for (e in seq_len(nrow(mcoords))){
        s <- mcoords[e,]
        #print (s)
        
        s1 <- data.frame(f1 = s$V1, f2 = s$n + block_len, f3 = nn, f4 = 1)
        s2 <- data.frame(f1 = s$V2, f2 = s$n + block_len, f3 = nn, f4 = 2)
        s3 <- data.frame(f1 = s$V4, f2 = s$n - block_len, f3 = nn, f4 = 3)
        s4 <- data.frame(f1 = s$V3, f2 = s$n - block_len, f3 = nn, f4 = 4)
        
        if (file.info(paste0(prefix, ".snps"))$size != 0){
            snps_t1 <- as.numeric(strsplit(strsplit(as.character(snps$V11[1]),":")[[1]][4],"-")[[1]][1])
            snps_t2 <- as.numeric(strsplit(strsplit(as.character(snps$V12[1]),":")[[1]][4],"-")[[1]][1])
            snps$st <- s$n + block_len
            snps$ed <- s$n - block_len
        }else{
            TRUE
        }
        
        if (nn == 1){
            sout <- rbind(s1,s2,s3,s4)
            if (file.info(paste0(prefix, ".snps"))$size != 0){
                snpst <- snps
            }else{
                TRUE
            }
        }else{
            sout <- rbind(sout, s1, s2, s3, s4)
            if (file.info(paste0(prefix, ".snps"))$size != 0){
                snpst <- rbind(snpst, snps)
            }else{
                TRUE
            }
        }
        
        nn = nn +1
        
    }
    
    #print (sout)
}


id <- unique(gfftotal[,c("gene_h","prefix")])

p <- ggplot() + geom_polygon(sout, mapping=aes(x = f1, y = f2, group  = f3 ), linetype = 1, fill = "lightblue")+ 
      #geom_point(s, mapping = aes(x = V1, y = 5), color = "red", size = 0.2) + 
      #geom_point(s, mapping = aes(x = V4, y = 0), color = "red", size = 0.2) + 
      #geom_segment(aes(x = V1, y = st, xend = V4, yend = ed), alpha = .8, color = "#F2735E", size = .5 , data = snpst) + 
      geom_segment(data=gfftotal[gfftotal$V6=="gene",],
                   mapping = aes(x = V7-V2, y = (exon_h1 + exon_h2)/2, xend = V8-V2, yend=(exon_h1 + exon_h2)/2), 
                   arrow = arrow(length=unit(0.15, "cm"), ends = gfftotal[gfftotal$V6=="gene",]$arrow)) + 
      geom_rect(aes(xmin = V7-V2, ymin = exon_h1, xmax = V8-V2, ymax = exon_h2), colour = NA, alpha = .4, size=.1 , data = gfftotal[gfftotal$V6=="exon",]) +
      #geom_vline(xintercept = c(54542107-54540000,54605188-54540000)) +
      theme_minimal() + coord_cartesian(xlim = c(0+zoomin,max(max(snpst$V7),max(snpst$V8))-zoomin)) + 
      scale_y_continuous(breaks = id$gene_h, labels = id$prefix) + ylab("") + xlab("")
#scale_x_continuous(limits =  c(0,max(max(snpst$V7),max(snpst$V8))), breaks = c(0,max(max(snpst$V7),max(snpst$V8))))# + ggtitle(prefix)

if (snpindel == TRUE){p <- p + geom_segment(aes(x = V1, y = st, xend = V4, yend = ed), alpha = .1, color = "#F2735E", size = .05 , data = snpst)}

pdf(paste0(prefix, ".alignment.pdf"), width = ww, height = hh)
print (p)
dev.off()


#print (p)
