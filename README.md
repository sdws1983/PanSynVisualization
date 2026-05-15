# PanSynVisualization

Pan-genome Synteny Map Visualization — extract, align, and visualize synteny between multiple genomes with optional phylogenetic tree integration.

## Quick start

```bash
Rscript make_synteny_vis.r --info example/test.info --out my_output --run
```

## Requirements

### External tools (must be on `$PATH`)

| Tool | Version | Purpose |
|------|---------|---------|
| [bedtools](https://bedtools.readthedocs.io/) | ≥2.30 | Extract sequences by region, intersect GFF annotations |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | ≥2.3 | Reverse-complement FASTA when strand is negative |
| [MUMmer4](https://mummer4.github.io/) (dnadiff) | ≥4.0 | Pairwise alignment (required for `--aligner mummer`) |
| [ParaFly](https://parafly.sourceforge.net/) | any | Parallelize dnadiff (required for `--aligner mummer`) |
| [mashtree](https://github.com/lskatz/mashtree) | any | *(optional)* Build phylogenetic tree (`--tree auto`) |

### R packages

| Package | Purpose |
|---------|---------|
| argparse | CLI argument parsing |
| ggplot2 | Synteny visualization |
| ape | Phylogenetic tree I/O |
| ggtree | Tree visualization (ggplot2-based) |
| patchwork | Combine tree + synteny plots |
| [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) | *(optional)* Pure-R local alignment (`--aligner biostrings`) |

## Usage

```
Rscript make_synteny_vis.r --info <file> --out <prefix> [options]
```

### Required arguments

| Flag | Description |
|------|-------------|
| `--info` | Genomic information file (tab-separated, see below) |
| `--out` | Output file prefix |

### Optional arguments

| Flag | Default | Description |
|------|---------|-------------|
| `--upstream` | 10000 | Upstream padding (bp) added before each region |
| `--downstream` | 10000 | Downstream padding (bp) added after each region |
| `--zoomin` | 0 | Zoom in from right side (bp) |
| `--highlight` | — | File with gene names to highlight (one per line) |
| `--outdir` | `.` | Output directory for the final PDF |
| `--width` | 6 | Figure width in inches |
| `--height` | 3 | Figure height in inches |
| `--snp` | TRUE | Show SNP/indel tracks on synteny blocks |
| `--run` | TRUE | Execute external commands; set `FALSE` for dry-run |
| `--cpu` | auto | CPU cores for ParaFly (default: detected cores − 1) |
| `--keep-tmp` | FALSE | Keep intermediate files for debugging |
| `--tree` | — | Phylogenetic tree: `auto` or path to `.nwk` file |
| `--tree-width` | 0.35 | Relative width of the tree panel |
| `--aligner` | mummer | Alignment engine: `mummer` or `biostrings` |

## Input format

The `--info` file is **tab-separated** with these columns:

| Column | Description |
|--------|-------------|
| `Prefix` | Identifier used in output filenames |
| `Genome` | Path to reference genome FASTA |
| `Region` | Genomic region in `Chr:start-end` format |
| `GFF` | Path to GFF3 annotation file |
| `Strand` | `TRUE` (forward) or `FALSE` (reverse complement) |

Example (`example/test.info`):

```
Prefix	Genome	Region	GFF	Strand
Npp.06A001160.1	~/ref/genome.fasta	Chr6A:3851740-3863740	~/ref/annot.gff3	TRUE
Npp.06B001620.1	~/ref/genome.fasta	Chr6B:4466259-4478259	~/ref/annot.gff3	TRUE
Npp.06C001520.1	~/ref/genome.fasta	Chr6C:5943903-5955903	~/ref/annot.gff3	TRUE
Npp.06D001220.1	~/ref/genome.fasta	Chr6D:4194613-4220850	~/ref/annot.gff3	TRUE
```

Each consecutive pair of rows produces one alignment plot (row 1 vs row 2, row 2 vs row 3, …).

## Pipeline stages

1. **Extract** — `bedtools getfasta` slices each region (± upstream/downstream) from the genome FASTA. If `Strand=FALSE`, `seqkit` reverse-complements the sequence. `bedtools intersect` extracts gene/exon annotations from the GFF3.

2. **Align** — Pairwise alignment for each consecutive pair. With `--aligner mummer` (default): `dnadiff` + `ParaFly` in parallel. With `--aligner biostrings`: R-based local alignment with iterative block extraction.

3. **Plot** — Alignment coordinates are rendered as synteny block polygons. Gene arrows and exon rectangles are drawn from GFF annotations. SNP/indel tracks are optionally overlaid. Output is a PDF.

## Alignment engines

Two alignment backends are available via `--aligner`:

### MUMmer (`--aligner mummer`, default)

Uses `dnadiff` (MUMmer4) and `ParaFly` for fast, parallel pairwise alignment. Finds discrete synteny blocks using maximal unique matches (MUMs) — ideal for cross-species comparisons and detecting inversions, duplications, and rearrangements.

### Biostrings (`--aligner biostrings`)

Pure-R implementation using the [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) package. Performs local pairwise alignment with iterative masking to extract synteny blocks at base-pair resolution. No external tools needed beyond bedtools and seqkit.

```bash
Rscript make_synteny_vis.r --info test.info --out demo --aligner biostrings --run
```

Characteristics:
- Base-pair-level block boundaries (every gap in the alignment is resolved)
- Works entirely in R — no MUMmer or ParaFly installation needed
- Slower than MUMmer for large sequences (~20–30s per 30kb pair)
- Can be combined with `--tree` mode
- SNP/indel tracks are not available (Biostrings does not produce `.snps` files)

## Phylogenetic tree mode

Add `--tree` to integrate a phylogenetic tree alongside the synteny plot:

```bash
# Build tree automatically from extracted sequences with mashtree
Rscript make_synteny_vis.r --info test.info --out demo --tree auto --run

# Use a pre-built Newick tree
Rscript make_synteny_vis.r --info test.info --out demo --tree my_tree.nwk --run
```

When a tree is provided:
- Genome rows are reordered to match the tree topology
- The tree is displayed on the **left** side of the plot
- Synteny blocks are shown on the right, y-axis aligned with tree tips
- `Prefix` names appear as y-axis labels

The tree tip labels must match the `Prefix` column in the info file. Extra rows not present in the tree are placed at the bottom.

## Features

- **Resume support** — Re-running the same command skips already-completed extraction, annotation, and alignment steps
- **Dry-run mode** — `--run FALSE` prints all commands without executing them
- **Dependency checking** — Verifies all external tools are on `PATH` before starting
- **Structured logging** — Timestamped `[INFO]`/`[WARN]`/`[ERROR]` messages
- **Input validation** — Checks column names, region format, and file existence before processing
- **Automatic cleanup** — Temporary files are removed on completion (use `--keep-tmp` to preserve)

## Output

The pipeline produces a PDF file: `<out>.alignment.pdf`

Example outputs are in the `example/` directory.
