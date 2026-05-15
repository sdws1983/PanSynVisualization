# PanSynVisualization — Pan-genome Synteny Map Visualization
# Entry point: Rscript make_synteny_vis.r --info <info_file> --out <prefix> [--run]
library(argparse)
library(ggplot2)
library(ape)
library(ggtree)
library(patchwork)

# ── Logging helpers ──────────────────────────────────────────────────────────

log_info  <- function(...) cat(paste0("[", Sys.time(), "] INFO  ", ..., "\n"))
log_warn  <- function(...) cat(paste0("[", Sys.time(), "] WARN  ", ..., "\n"))
log_error <- function(...) cat(paste0("[", Sys.time(), "] ERROR ", ..., "\n"))

# ── Safe system call ─────────────────────────────────────────────────────────

system_or_stop <- function(cmd, msg = cmd, ignore_stderr = FALSE) {
  log_info("Running: ", cmd)
  stderr_redir <- if (ignore_stderr) " 2>/dev/null" else ""
  rc <- system(paste0(cmd, stderr_redir))
  if (rc != 0) stop(msg, " (exit code ", rc, ")", call. = FALSE)
  invisible(rc)
}

# ── Dependency checker ──────────────────────────────────────────────────────

check_dependencies <- function() {
  required <- c(bedtools = "bedtools", seqkit = "seqkit", dnadiff = "dnadiff", ParaFly = "ParaFly")
  version_flags <- c(
    bedtools = "--version",
    seqkit   = "version",
    dnadiff  = "--version",
    ParaFly  = "--version"
  )
  for (tool in names(required)) {
    exe <- required[tool]
    if (Sys.which(exe) == "")
      stop("Required tool '", exe, "' not found on PATH. Please install it.", call. = FALSE)

    ver <- tryCatch({
      out <- system(paste(exe, version_flags[tool], "2>&1 | head -1"), intern = TRUE)
      trimws(paste(out, collapse = " "))
    }, error = function(e) "available")
    log_info("Found ", exe, ": ", ver)
  }
}

# ── Tree helpers ─────────────────────────────────────────────────────────────

build_tree_mashtree <- function(fa_files, tmpdir, cpu) {
  if (Sys.which("mashtree") == "")
    stop("mashtree not found on PATH. Install it: conda install -c bioconda mashtree", call. = FALSE)

  log_info("Building Mash tree from ", length(fa_files), " sequences...")
  fa_list <- paste(fa_files, collapse = " ")
  tree_path <- file.path(tmpdir, "pansyn_tree.nwk")
  cmd <- paste0("mashtree --numcpus ", cpu, " ", fa_list, " > ", tree_path)
  system_or_stop(cmd, "mashtree failed")
  log_info("Tree saved to: ", tree_path)
  return(tree_path)
}

read_and_validate_tree <- function(tree_path, prefixes) {
  if (!file.exists(tree_path))
    stop("Tree file not found: ", tree_path, call. = FALSE)

  tree <- ape::read.tree(tree_path)
  tree <- ape::ladderize(tree)

  # Tree tip labels must be a subset of info Prefixes
  missing <- setdiff(tree$tip.label, prefixes)
  extra   <- setdiff(prefixes, tree$tip.label)
  if (length(missing) > 0)
    stop("Tree tip(s) not found in info file: ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(extra) > 0)
    log_warn("Info row(s) not in tree (will be placed at end): ", paste(extra, collapse = ", "))

  # ggtree places tips at y=1 (bottom) to y=N (top)
  # We want tips in tip.label order from bottom to top
  log_info("Tree loaded: ", length(tree$tip.label), " tips")
  return(tree)
}

reorder_info_by_tree <- function(info, tree) {
  # Reorder info rows to match tree tip order (bottom-to-top)
  # Tips not in the tree go at the end (bottom)
  tip_order <- tree$tip.label
  in_tree    <- info$Prefix %in% tip_order
  info_tree  <- info[in_tree, ]
  info_other <- info[!in_tree, ]

  # Sort info_tree rows by tree tip order
  info_tree <- info_tree[match(tip_order, info_tree$Prefix), ]
  info_tree <- info_tree[!is.na(info_tree$Prefix), ]

  info_new <- rbind(info_other, info_tree)  # extras at bottom
  log_info("Genomes reordered by tree: ", paste(info_new$Prefix, collapse = ", "))
  return(info_new)
}

# ── Input validation ─────────────────────────────────────────────────────────

validate_info <- function(info) {
  required_cols <- c("Prefix", "Genome", "Region", "GFF", "Strand")
  missing_cols <- setdiff(required_cols, colnames(info))
  if (length(missing_cols) > 0)
    stop("Missing required columns in info file: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  if (nrow(info) < 2)
    stop("Info file must contain at least 2 rows (need a pair for alignment)", call. = FALSE)

  for (i in seq_len(nrow(info))) {
    region <- info[i,]$Region
    if (!grepl("^[^:]+:[0-9]+-[0-9]+$", region))
      stop("Row ", i, ": Region must be 'Chr:start-end', got: ", region, call. = FALSE)
    if (!file.exists(info[i,]$Genome))
      stop("Row ", i, ": Genome file not found: ", info[i,]$Genome, call. = FALSE)
    if (!file.exists(info[i,]$GFF))
      stop("Row ", i, ": GFF file not found: ", info[i,]$GFF, call. = FALSE)
  }
  log_info("Info file validated: ", nrow(info), " rows, ", nrow(info) - 1, " pairwise alignment(s)")
}

# ── Extract sequence ─────────────────────────────────────────────────────────

extract_sequence <- function(chro, start, end, genome, fa_path, strand, run) {
  cmd <- sprintf("printf '%s\t%s\t%s\n' %s %s %s | bedtools getfasta -bed - -fi %s -fo %s",
                 chro, start, end, chro, start, end, genome, fa_path)

  if (!run) {
    log_info("[DRY-RUN] Would run: ", cmd)
    return(invisible())
  }

  if (file.exists(fa_path) && file.info(fa_path)$size > 0) {
    log_info("Skipping (exists): ", fa_path)
    return(invisible())
  }

  system_or_stop(cmd, paste0("bedtools getfasta failed for ", fa_path))

  if (identical(strand, FALSE) || identical(strand, "FALSE")) {
    log_info("Reverse-complementing: ", fa_path)
    cmd_rc <- paste0("seqkit seq -t dna -rp ", fa_path, " > ", fa_path, ".rc && mv ", fa_path, ".rc ", fa_path)
    system_or_stop(cmd_rc, paste0("seqkit reverse-complement failed for ", fa_path))
  }
}

# ── Extract annotations ──────────────────────────────────────────────────────

extract_annotations <- function(chro, start, end, gff, gff3_path, strand, nn, exon_h, run) {
  cmd <- sprintf("printf '%s\t%s\t%s\n' %s %s %s | bedtools intersect -a - -b %s -wo > %s",
                 chro, start, end, chro, start, end, gff, gff3_path)

  if (!run) {
    log_info("[DRY-RUN] Would run: ", cmd)
    return(data.frame())
  }

  if (file.exists(gff3_path) && file.info(gff3_path)$size > 0) {
    log_info("Skipping (exists): ", gff3_path)
  } else {
    system_or_stop(cmd, paste0("bedtools intersect failed for ", gff3_path))
  }

  gf <- data.frame()
  if (file.exists(gff3_path) && file.info(gff3_path)$size > 0) {
    gf <- read.table(gff3_path, header = FALSE, sep = "\t")
    gf <- gf[(gf$V6 == "gene") | (gf$V6 == "exon"), ]
    gf$exon_h1 <- nn - exon_h
    gf$exon_h2 <- nn + exon_h
    gf$gene_h  <- nn
    if (identical(strand, FALSE) || identical(strand, "FALSE")) {
      gf$V8 <- gf$V3 - gf$V8 + gf$V2
      gf$V7 <- gf$V3 - gf$V7 + gf$V2
    }
  } else {
    log_warn("No gene/exon annotations found for: ", gff3_path)
  }
  return(gf)
}

# ── Run alignments ───────────────────────────────────────────────────────────

run_alignments <- function(out, tmpdir, command_file, cpu, run) {
  # Build pending commands (skip if .mcoords already exists)
  pending <- character(0)
  for (i in seq_len(nrow(out))) {
    l <- out[i, ]
    prefix      <- paste0(l$prefix1, "_", l$prefix2)
    mcoords_path <- file.path(tmpdir, paste0(prefix, ".mcoords"))
    if (!run) {
      fa1_abs <- file.path(tmpdir, basename(l$fa1))
      fa2_abs <- file.path(tmpdir, basename(l$fa2))
      pending <- c(pending, paste0("dnadiff -p ", file.path(tmpdir, prefix), " ", fa1_abs, " ", fa2_abs))
    } else if (file.exists(mcoords_path) && file.info(mcoords_path)$size > 0) {
      log_info("Skipping alignment (exists): ", basename(mcoords_path))
    } else {
      fa1_abs <- file.path(tmpdir, basename(l$fa1))
      fa2_abs <- file.path(tmpdir, basename(l$fa2))
      pending <- c(pending, paste0("dnadiff -p ", file.path(tmpdir, prefix), " ", fa1_abs, " ", fa2_abs))
    }
  }

  writeLines(pending, command_file)

  if (!run) {
    log_info("Dry-run: wrote ", length(pending), " alignment commands to ", command_file)
    return()
  }
  if (length(pending) == 0) {
    log_info("All alignments already complete")
    return()
  }

  log_info("Running ", length(pending), " alignment(s) with ParaFly (CPU=", cpu, ")")
  cmd <- paste0("ParaFly -c ", command_file, " -CPU ", cpu)
  system_or_stop(cmd, "ParaFly alignment step failed")

  # Verify all expected outputs exist
  for (i in seq_len(nrow(out))) {
    l <- out[i, ]
    mcoords_path <- file.path(tmpdir, paste0(l$prefix1, "_", l$prefix2, ".mcoords"))
    if (!file.exists(mcoords_path))
      stop("Missing expected output: ", mcoords_path, call. = FALSE)
  }
  log_info("All ", nrow(out), " alignments completed successfully")
}

# ── Parse alignment results ──────────────────────────────────────────────────

parse_alignments <- function(out, tmpdir, snpindel) {
  mcoords_list <- list()
  snps_list    <- list()

  for (i in seq_len(nrow(out))) {
    l      <- out[i, ]
    prefix <- paste0(l$prefix1, "_", l$prefix2)
    mcoords_path <- file.path(tmpdir, paste0(prefix, ".mcoords"))

    if (!file.exists(mcoords_path))
      stop("Alignment output missing: ", mcoords_path, call. = FALSE)

    mcoords <- read.table(mcoords_path, header = FALSE)
    mcoords$prefix <- prefix
    mcoords_list[[i]] <- mcoords

    if (snpindel) {
      snps_path <- file.path(tmpdir, paste0(prefix, ".snps"))
      if (file.exists(snps_path) && file.info(snps_path)$size > 0) {
        snps <- read.table(snps_path, header = FALSE)
        snps$prefix <- prefix
        snps_list[[i]] <- snps
      } else {
        snps_list[[i]] <- data.frame()
      }
    } else {
      snps_list[[i]] <- data.frame()
    }
  }
  list(mcoords = mcoords_list, snps = snps_list)
}

# ── Build plot data ──────────────────────────────────────────────────────────

build_synteny_polygons <- function(mcoords_list, snps_list, nn_start = 1.5, block_len = 0.40) {
  sout  <- data.frame()
  snpst <- data.frame()
  nn <- 1

  for (i in seq_along(mcoords_list)) {
    mcoords <- mcoords_list[[i]]
    snps    <- snps_list[[i]]
    # pair_n sits halfway between adjacent genome rows:
    # gene_h[j] = nn_start - (j-1), so pair between j and j+1 is at nn_start - j + 0.5
    pair_n  <- nn_start - i + 0.5

    for (e in seq_len(nrow(mcoords))) {
      s <- mcoords[e, ]
      s1 <- data.frame(f1 = s$V1, f2 = pair_n + block_len, f3 = nn)
      s2 <- data.frame(f1 = s$V2, f2 = pair_n + block_len, f3 = nn)
      s3 <- data.frame(f1 = s$V4, f2 = pair_n - block_len, f3 = nn)
      s4 <- data.frame(f1 = s$V3, f2 = pair_n - block_len, f3 = nn)

      if (nn == 1) {
        sout <- rbind(s1, s2, s3, s4)
      } else {
        sout <- rbind(sout, s1, s2, s3, s4)
      }

      if (nrow(snps) > 0) {
        snps$st <- pair_n + block_len
        snps$ed <- pair_n - block_len
        if (nn == 1) {
          snpst <- snps
        } else {
          snpst <- rbind(snpst, snps)
        }
      }
      nn <- nn + 1
    }
  }
  list(sout = sout, snpst = snpst)
}

# ── Render plot ──────────────────────────────────────────────────────────────

render_plot <- function(sout, gfftotal, snpst, snpindel, zoomin, prefixo, outdir, ww, hh) {
  id <- unique(gfftotal[, c("gene_h", "prefix")])

  p <- ggplot() +
    geom_polygon(sout, mapping = aes(x = f1, y = f2, group = f3),
                 linetype = 1, fill = "lightblue") +
    geom_segment(
      data    = gfftotal[gfftotal$V6 == "gene", ],
      mapping = aes(x = V7 - V2, y = (exon_h1 + exon_h2) / 2,
                    xend = V8 - V2, yend = (exon_h1 + exon_h2) / 2,
                    color = highlight),
      arrow   = arrow(length = unit(0.15, "cm"),
                      ends = gfftotal[gfftotal$V6 == "gene", ]$arrow)
    ) +
    geom_rect(
      data    = gfftotal[gfftotal$V6 == "exon", ],
      mapping = aes(xmin = V7 - V2, ymin = exon_h1,
                    xmax = V8 - V2, ymax = exon_h2),
      colour  = NA, alpha = 0.4, size = 0.1
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, max(0, max(gfftotal$V8 - gfftotal$V2) - zoomin))) +
    scale_y_continuous(breaks = id$gene_h, labels = id$prefix) +
    ylab("") + xlab("")

  if (snpindel && nrow(snpst) > 0) {
    p <- p + geom_segment(
      aes(x = V1, y = st, xend = V4, yend = ed),
      alpha = 0.1, color = "#F2735E", size = 0.05, data = snpst
    )
  }

  pdf_path <- file.path(outdir, paste0(prefixo, ".alignment.pdf"))
  pdf(pdf_path, width = ww, height = hh)
  print(p)
  dev.off()
  log_info("Plot saved to: ", pdf_path)
}

# ── Render plot with tree ─────────────────────────────────────────────────────

render_plot_with_tree <- function(sout, gfftotal, snpst, snpindel, zoomin,
                                   prefixo, outdir, ww, hh, tree, tree_width) {
  id <- unique(gfftotal[, c("gene_h", "prefix")])

  # Synteny plot
  p_syn <- ggplot() +
    geom_polygon(sout, mapping = aes(x = f1, y = f2, group = f3),
                 linetype = 1, fill = "lightblue") +
    geom_segment(
      data    = gfftotal[gfftotal$V6 == "gene", ],
      mapping = aes(x = V7 - V2, y = (exon_h1 + exon_h2) / 2,
                    xend = V8 - V2, yend = (exon_h1 + exon_h2) / 2,
                    color = highlight),
      arrow   = arrow(length = unit(0.15, "cm"),
                      ends = gfftotal[gfftotal$V6 == "gene", ]$arrow)
    ) +
    geom_rect(
      data    = gfftotal[gfftotal$V6 == "exon", ],
      mapping = aes(xmin = V7 - V2, ymin = exon_h1,
                    xmax = V8 - V2, ymax = exon_h2),
      colour  = NA, alpha = 0.4, size = 0.1
    ) +
    theme_minimal() +
    scale_y_continuous(breaks = id$gene_h, labels = id$prefix) +
    ylab("") + xlab("")

  if (snpindel && nrow(snpst) > 0) {
    p_syn <- p_syn + geom_segment(
      aes(x = V1, y = st, xend = V4, yend = ed),
      alpha = 0.1, color = "#F2735E", size = 0.05, data = snpst
    )
  }

  # Tree plot — tips point right toward synteny, labels match Prefix
  n_genomes <- length(tree$tip.label)
  p_tree <- ggtree(tree, ladderize = TRUE, size = 0.3) +
    theme_tree2() +
    ylim(0.5, n_genomes + 0.5) +
    theme(plot.margin = margin(0, 0, 0, 5))

  # Set shared x/y limits on synteny (use coord_cartesian to avoid overriding scale)
  x_max <- max(0, max(gfftotal$V8 - gfftotal$V2) - zoomin)
  p_syn <- p_syn + coord_cartesian(xlim = c(0, x_max), ylim = c(0.5, n_genomes + 0.5))

  # Combine: tree (left) + synteny (right)
  p_combined <- wrap_plots(
    p_tree,
    p_syn + theme(plot.margin = margin(5, 5, 5, 0)),
    nrow = 1,
    widths = c(tree_width, 1 - tree_width)
  )

  pdf_path <- file.path(outdir, paste0(prefixo, ".alignment.pdf"))
  pdf(pdf_path, width = ww, height = hh)
  print(p_combined)
  dev.off()
  log_info("Plot with tree saved to: ", pdf_path)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Main pipeline
# ═══════════════════════════════════════════════════════════════════════════════

main <- function() {
  # ── Parse arguments ──────────────────────────────────────────────────────
  parser <- ArgumentParser(description = "Pangenome Synteny visualization")
  parser$add_argument("--info",       help = "Genomic information file", required = TRUE)
  parser$add_argument("--upstream",   type = "numeric", help = "Upstream padding (default=10000)",  default = 10000)
  parser$add_argument("--downstream", type = "numeric", help = "Downstream padding (default=10000)", default = 10000)
  parser$add_argument("--zoomin",     type = "numeric", help = "Zoom in from right (default=0)",     default = 0)
  parser$add_argument("--highlight",  help = "File with gene names to highlight", default = FALSE)
  parser$add_argument("--out",        help = "Output file prefix", required = TRUE)
  parser$add_argument("--outdir",     help = "Output directory (default: current)", default = ".")
  parser$add_argument("--width",      type = "numeric", help = "Figure width in inches (default=6)",  default = 6)
  parser$add_argument("--height",     type = "numeric", help = "Figure height in inches (default=3)", default = 3)
  parser$add_argument("--snp",        help = "Show SNP/indel tracks",   default = TRUE)
  parser$add_argument("--run",        help = "Actually execute external commands", default = TRUE)
  parser$add_argument("--cpu",        type = "integer", help = "CPU cores for ParaFly", default = NULL)
  parser$add_argument("--keep-tmp",   help = "Keep temporary files", default = FALSE)
  parser$add_argument("--tree",       help = "Phylogenetic tree: 'auto' to build with mashtree, or path to .nwk file", default = FALSE)
  parser$add_argument("--tree-width", type = "numeric", help = "Relative width of tree panel (default=0.35)", default = 0.35)

  args <- parser$parse_args()

  run       <- as.logical(args$run)
  snpindel  <- as.logical(args$snp)
  keep_tmp  <- as.logical(args$keep_tmp)
  cpu       <- if (is.null(args$cpu)) max(1, parallel::detectCores() - 1) else args$cpu
  upstream  <- args$upstream
  downstream <- args$downstream
  zoomin    <- args$zoomin
  prefixo   <- args$out
  outdir    <- args$outdir
  ww        <- args$width
  hh        <- args$height
  highlightgene <- args$highlight
  tree_mode  <- args$tree
  tree_width <- args$tree_width

  # ── Setup ────────────────────────────────────────────────────────────────
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  tmpdir <- file.path(outdir, paste0(".pansyn_", prefixo))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

  on_exit_cleanup <- function() {
    if (!keep_tmp) {
      unlink(tmpdir, recursive = TRUE)
      log_info("Cleaned up temporary directory")
    } else {
      log_info("Temporary files kept at: ", tmpdir)
    }
  }
  on.exit(on_exit_cleanup(), add = TRUE)

  log_info("PanSynVisualization starting")
  log_info("Temporary directory: ", tmpdir)
  log_info("Output directory: ", outdir)
  log_info("CPU cores: ", cpu)

  check_dependencies()

  # ── Read and validate input ────────────────────────────────────────────
  info <- read.table(args$info, header = TRUE)
  validate_info(info)

  highlight_genes <- NULL
  if (!identical(highlightgene, FALSE) && file.exists(highlightgene)) {
    highlight_genes <- read.table(highlightgene, header = FALSE)$V1
    log_info("Loaded ", length(highlight_genes), " highlight genes")
  }

  # ── Constants ──────────────────────────────────────────────────────────
  block_len  <- 0.40
  exon_h     <- 0.06
  using_tree <- !identical(tree_mode, FALSE)
  tree       <- NULL

  # ── Pre-extraction for tree auto mode ────────────────────────────────
  if (identical(tree_mode, "auto")) {
    log_info("── Pre-extraction: extracting sequences for tree building ──")
    fa_files <- character(nrow(info))
    for (i in seq_len(nrow(info))) {
      row   <- info[i, ]
      region <- row$Region
      t     <- strsplit(region, ":")[[1]]
      chro  <- t[1]
      se    <- as.numeric(strsplit(t[2], "-")[[1]])
      start <- se[1] - upstream
      end   <- se[2] + downstream
      fa_filename <- paste0(row$Prefix, "_", chro, "_", start, "_", end, ".fasta")
      fa_path <- file.path(tmpdir, fa_filename)
      extract_sequence(chro, start, end, row$Genome, fa_path, row$Strand, run)
      fa_files[i] <- fa_path
    }
    tree_path <- build_tree_mashtree(fa_files[file.exists(fa_files)], tmpdir, cpu)
    tree_mode <- tree_path  # now use this tree file for the rest
  }

  # ── Load and validate tree ────────────────────────────────────────────
  if (using_tree) {
    log_info("── Loading phylogenetic tree ──")
    tree <- read_and_validate_tree(tree_mode, info$Prefix)
    info <- reorder_info_by_tree(info, tree)
    # Reverse so first row = top tip (gene_h = N), last row = bottom tip (gene_h = 1)
    info <- info[nrow(info):1, ]
    rownames(info) <- NULL
  }

  # ── Stage 1: Extract sequences and annotations ────────────────────────
  log_info("── Stage 1: Extracting sequences and annotations ──")

  n         <- 1
  gfftotal  <- data.frame()
  out       <- NULL
  l_fa <- l_chro <- l_start <- l_end <- l_pre <- NULL
  # In tree mode, nn starts at N (top genome); otherwise starts at 1.5
  nn <- if (using_tree) nrow(info) else 1.5

  for (i in seq_len(nrow(info))) {
    row    <- info[i, ]
    region <- row$Region
    t      <- strsplit(region, ":")[[1]]
    chro   <- t[1]
    se     <- as.numeric(strsplit(t[2], "-")[[1]])
    start  <- se[1] - upstream
    end    <- se[2] + downstream

    fa_filename <- paste0(row$Prefix, "_", chro, "_", start, "_", end, ".fasta")
    fa_path     <- file.path(tmpdir, fa_filename)

    log_info(sprintf("Row %d/%d: %s (%s:%d-%d)", i, nrow(info), row$Prefix, chro, start, end))

    # Extract sequence (resume-safe: skips if fa_path exists from pre-extraction)
    extract_sequence(chro, start, end, row$Genome, fa_path, row$Strand, run)

    # Extract annotations
    gff3_filename <- paste0(row$Prefix, ".gff3")
    gff3_path     <- file.path(tmpdir, gff3_filename)
    gf <- extract_annotations(chro, start, end, row$GFF, gff3_path, row$Strand, nn, exon_h, run)
    if (nrow(gf) > 0) gf$prefix <- row$Prefix

    nn <- nn - 1

    # Build pairwise output table
    if (n == 1) {
      l_fa <- fa_path; l_chro <- chro; l_start <- start; l_end <- end; l_pre <- row$Prefix
      n <- 0
      gfftotal <- gf
    } else if (n == 0) {
      out <- data.frame(
        prefix1 = l_pre, prefix2 = row$Prefix,
        fa1 = l_fa, fa2 = fa_path,
        chr1 = l_chro, chr2 = chro,
        start1 = l_start, start2 = start,
        end1 = l_end, end2 = end,
        stringsAsFactors = FALSE
      )
      l_fa <- fa_path; l_chro <- chro; l_start <- start; l_end <- end; l_pre <- row$Prefix
      n <- -1
      gfftotal <- rbind(gfftotal, gf)
    } else {
      out <- rbind(out, data.frame(
        prefix1 = l_pre, prefix2 = row$Prefix,
        fa1 = l_fa, fa2 = fa_path,
        chr1 = l_chro, chr2 = chro,
        start1 = l_start, start2 = start,
        end1 = l_end, end2 = end,
        stringsAsFactors = FALSE
      ))
      l_fa <- fa_path; l_chro <- chro; l_start <- start; l_end <- end; l_pre <- row$Prefix
      gfftotal <- rbind(gfftotal, gf)
    }
  }

  log_info(sprintf("Extracted %d regions, %d alignment pairs", nrow(info), nrow(out)))

  # ── Stage 2: Align ────────────────────────────────────────────────────
  log_info("── Stage 2: Pairwise alignment ──")

  command_file <- file.path(tmpdir, "align_commands.txt")
  run_alignments(out, tmpdir, command_file, cpu, run)

  # ── Stage 3: Build gene arrows and highlights ────────────────────────
  log_info("── Stage 3: Building plot data ──")

  if (nrow(gfftotal) > 0) {
    gfftotal$arrow <- ""
    gfftotal[gfftotal$V10 == "+", ]$arrow <- "last"
    gfftotal[gfftotal$V10 == "-", ]$arrow <- "first"
  }

  if (nrow(gfftotal) > 0) {
    gfftotal$highlight <- "no"
    if (!is.null(highlight_genes)) {
      gfftotal[gfftotal$V12 %in% highlight_genes, ]$highlight <- "highlight"
    }
  }

  # ── Stage 4: Parse alignment results and build polygons ───────────────

  if (!run) {
    log_info("Dry-run complete — stopping before parse/plot. Use --run TRUE to execute fully.")
    return(invisible(NULL))
  }

  log_info("── Stage 4: Parsing alignment results ──")

  results <- parse_alignments(out, tmpdir, snpindel)
  nn_start <- if (using_tree) nrow(info) else 1.5
  poly    <- build_synteny_polygons(results$mcoords, results$snps, nn_start)

  # ── Stage 5: Render plot ──────────────────────────────────────────────
  log_info("── Stage 5: Rendering plot ──")

  if (using_tree && !is.null(tree)) {
    render_plot_with_tree(poly$sout, gfftotal, poly$snpst, snpindel, zoomin,
                          prefixo, outdir, ww, hh, tree, tree_width)
  } else {
    render_plot(poly$sout, gfftotal, poly$snpst, snpindel, zoomin, prefixo, outdir, ww, hh)
  }

  log_info("PanSynVisualization finished successfully")
}

main()
