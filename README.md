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
| [MUMmer4](https://mummer4.github.io/) (dnadiff) | ≥4.0 | Pairwise genome alignment |
| [ParaFly](https://parafly.sourceforge.net/) | any | Parallelize dnadiff across pairs |
| [mashtree](https://github.com/lskatz/mashtree) | any | *(optional)* Build phylogenetic tree from sequences (`--tree auto`) |

### R packages

| Package | Purpose |
|---------|---------|
| argparse | CLI argument parsing |
| ggplot2 | Synteny visualization |
| ape | Phylogenetic tree I/O |
| ggtree | Tree visualization (ggplot2-based) |
| patchwork | Combine tree + synteny plots |

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

2. **Align** — `dnadiff` performs pairwise alignment for each consecutive pair. `ParaFly` runs alignments in parallel.

3. **Plot** — MUMmer coordinates are rendered as synteny block polygons. Gene arrows and exon rectangles are drawn from GFF annotations. SNP/indel tracks are optionally overlaid. Output is a PDF.

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
