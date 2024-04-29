# PanSynVisualization
Pan-genome Synteny Map Visualization

## Requirement

ggplot2

[Seqkit](https://bioinf.shenwei.me/seqkit/)

[Mummer4](https://mummer4.github.io/)

[Bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

[ParaFly](https://parafly.sourceforge.net/)

## Useage

```
$ Rscript make_synteny_vis.r -h
usage: make_synteny_vis.r [-h] --info INFO [--upstream UPSTREAM]
                          [--downstream DOWNSTREAM] [--zoomin ZOOMIN] --out
                          OUT [--width WIDTH] [--height HEIGHT] [--snp]
                          [--run]

Pangenome Synteny visualization

optional arguments:
  -h, --help            show this help message and exit
  --info INFO           Genomic information file
  --upstream UPSTREAM   Upstream (default=10000)
  --downstream DOWNSTREAM
                        Downstream (default=10000)
  --zoomin ZOOMIN       Zoomin (default=0)
  --highlight HIGHLIGHT
                        Highlight genes
  --out OUT             output file prefix
  --tmp TMP             tmp file
  --width WIDTH         Figure width (default=6)
  --height HEIGHT       Figure height (default=3)
  --snp SNP             Flag to indicate if the snp/indels should be showed
  --run RUN             Flag to indicate if the mummer should run
```

## Input

[Info file format](example/test.info):

```
Prefix    Genome    Region    GFF       Strand
Npp.06A001160.1       ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.fasta   Chr6A:3851740-3863740   ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.gff3  TRUE
Npp.06B001620.1       ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.fasta   Chr6B:4466259-4478259   ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.gff3  TRUE
Npp.06C001520.1       ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.fasta   Chr6C:5943903-5955903 ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.gff3    TRUE
Npp.06D001220.1       ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.fasta   Chr6D:4194613-4220850 ~/reference/Np-X/Np-X.3ddna.Chr.v20210804.gff3    TRUE
```

## Output

[pdf file1](example/Npp.06C001520.1_Npp.06D001220.1.alignment.pdf)
[pdf file2](example/Srufi2_Chr1_11853340_11893344_Srufi_Chr01_12140000_12180000.alignment.pdf)