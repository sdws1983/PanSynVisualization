# PanSynVisualization
Pan-genome Synteny Map Visualization


------------------
### 1. make multiple sequences alignments
 - Environment and dependancies:
 
      Python3+
  
      [MUSCLE](http://www.drive5.com/muscle)
  
      [Bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
  
 - Use [example files](example) for script `plot_syn.py`:
 
 ```
 python3 plot_syn.py \
    -i info -f fasta/ -g gff3/ -u 5000 -d 3000 -o hsp101_5k_3k -p hsp101_5k_3k
 ```
 
 - Details:
 ```
 python3 plot_syn.py --help
usage: plot_syn.py [-h] [-i INPUT] [-o OUTPUT] [-p PREFIX] [-f FASTA] [-g GFF] [-u UP] [-d DOWN] [-a {True,False}]

----

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input gene info
  -o OUTPUT, --output OUTPUT
                        output folder
  -p PREFIX, --prefix PREFIX
                        output prefix
  -f FASTA, --fasta FASTA
                        fasta folder
  -g GFF, --gff GFF     gff folder
  -u UP, --up UP        gene upstream (default=10000)
  -d DOWN, --down DOWN  gene downstream (default=5000)
  -a {True,False}, --all {True,False}
                        run one step (default=True)
 ```
 - outputs:
 
    [pdf file](example/hsp101_5k_3k.pdf)
 
    [muscle alignment results](example/hsp101_5k_3k.2.fa)
 
    [muscle output tree](example/hsp101_5k_3k.2.nwk)


------------------
### 2. Synteny Map Visualization
 - If you turn on the “one step” option when running `plot_syn.py`, then there is no need to run `plot_syn.r` here again
 - Environment and dependancies:
 
      R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
  
 - Use [example files](example) for script `plot_syn.r`:
 
 ```
 Rscript plot_syn.r \
    --nwk hsp101_5k_3k.2.nwk --gffdir hsp101_5k_3k/ --up 5000 --down 3000 \
    hsp101_5k_3k.2.fa hsp101_5k_3k
 ```
 
 - Details:
```
$ Rscript plot_syn.r --help
usage: plot_syn.r [-h] [--nwk NWK] [--gffdir GFFDIR] [--up UP] [--down DOWN]
                  fi fo

Syntenic visualization using provided msa file

positional arguments:
  fi               msa file
  fo               output prefix

optional arguments:
  -h, --help       show this help message and exit
  --nwk NWK        nwk tree
  --gffdir GFFDIR  gff file dir
  --up UP          upstream
  --down DOWN      downstream
```
 
------------------
### 0. Pre-processing
 - Prepare reference genome folder
 ```
 $ ls fasta/
Zm-A188-REFERENCE-KSU-1.0.fa        Zm-CML247-REFERENCE-NAM-1.0.fa      Zm-CML52-REFERENCE-NAM-1.0.fa.fai  Zm-HP301-REFERENCE-NAM-1.0.fa.fai  Zm-Ky21-REFERENCE-NAM-1.0.fa       Zm-Mo18W-REFERENCE-NAM-1.0.fa.fai  Zm-Oh7B-REFERENCE-NAM-1.0.fa               Zm-Tx303-REFERENCE-NAM-1.0.fa
Zm-B73-REFERENCE-NAM-5.0.fa         Zm-CML247-REFERENCE-NAM-1.0.fa.fai  Zm-CML69-REFERENCE-NAM-1.0.fa      Zm-Ia453-REFERENCE-FL-1.0.fa       Zm-Ky21-REFERENCE-NAM-1.0.fa.fai   Zm-Ms71-REFERENCE-NAM-1.0.fa       Zm-Oh7B-REFERENCE-NAM-1.0.fa.fai           Zm-Tx303-REFERENCE-NAM-1.0.fa.fai
Zm-B73-REFERENCE-NAM-5.0.fa.fai     Zm-CML277-REFERENCE-NAM-1.0.fa      Zm-CML69-REFERENCE-NAM-1.0.fa.fai  Zm-Il14H-REFERENCE-NAM-1.0.fa      Zm-M162W-REFERENCE-NAM-1.0.fa      Zm-Ms71-REFERENCE-NAM-1.0.fa.fai   Zm-P39-REFERENCE-NAM-1.0.fa                Zm-Tzi8-REFERENCE-NAM-1.0.fa
Zm-B97-REFERENCE-NAM-1.0.fa         Zm-CML277-REFERENCE-NAM-1.0.fa.fai  Zm-DK105-REFERENCE-TUM-1.0.fa      Zm-Il14H-REFERENCE-NAM-1.0.fa.fai  Zm-M162W-REFERENCE-NAM-1.0.fa.fai  Zm-NC350-REFERENCE-NAM-1.0.fa      Zm-P39-REFERENCE-NAM-1.0.fa.fai            Zm-Tzi8-REFERENCE-NAM-1.0.fa.fai
Zm-B97-REFERENCE-NAM-1.0.fa.fai     Zm-CML322-REFERENCE-NAM-1.0.fa      Zm-EP1-REFERENCE-TUM-1.0.fa        Zm-K0326Y-REFERENCE-SIPPE-1.0.fa   Zm-M37W-REFERENCE-NAM-1.0.fa       Zm-NC350-REFERENCE-NAM-1.0.fa.fai  Zm-PE0075-REFERENCE-TUM-1.0.fa             Zm-W22-REFERENCE-NRGENE-2.0.fa
Zm-CML103-REFERENCE-NAM-1.0.fa      Zm-CML322-REFERENCE-NAM-1.0.fa.fai  Zm-EP1-REFERENCE-TUM-1.0.fa.fai    Zm-Ki11-REFERENCE-NAM-1.0.fa       Zm-M37W-REFERENCE-NAM-1.0.fa.fai   Zm-NC358-REFERENCE-NAM-1.0.fa      Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.fa      Zm-W22-REFERENCE-NRGENE-2.0.fa.fai
Zm-CML103-REFERENCE-NAM-1.0.fa.fai  Zm-CML333-REFERENCE-NAM-1.0.fa      Zm-F7-REFERENCE-TUM-1.0.fa         Zm-Ki11-REFERENCE-NAM-1.0.fa.fai   Zm-Mo17-REFERENCE-CAU-1.0.fa       Zm-NC358-REFERENCE-NAM-1.0.fa.fai  Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.fa.fai  Zx-mexicana-REFERENCE-YAN-1.0.fa
Zm-CML228-REFERENCE-NAM-1.0.fa      Zm-CML333-REFERENCE-NAM-1.0.fa.fai  Zm-F7-REFERENCE-TUM-1.0.fa.fai     Zm-Ki3-REFERENCE-NAM-1.0.fa        Zm-Mo17-REFERENCE-CAU-1.0.fa.fai   Zm-Oh43-REFERENCE-NAM-1.0.fa       Zm-SK-REFERENCE-YAN-1.0.fa                 Zx-mexicana-REFERENCE-YAN-1.0.fa.fai
Zm-CML228-REFERENCE-NAM-1.0.fa.fai  Zm-CML52-REFERENCE-NAM-1.0.fa       Zm-HP301-REFERENCE-NAM-1.0.fa      Zm-Ki3-REFERENCE-NAM-1.0.fa.fai    Zm-Mo18W-REFERENCE-NAM-1.0.fa      Zm-Oh43-REFERENCE-NAM-1.0.fa.fai   Zm-SK-REFERENCE-YAN-1.0.fa.fai             Zx-PI566673-REFERENCE-YAN-1.0.fa
 ```
 - Prepare genome annotation gff3 file folder
 ```
 ls gff3/
Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz     Zm-CML333-REFERENCE-NAM-1.0_Zm00026ab.1.gff3.gz  Zm-Il14H-REFERENCE-NAM-1.0_Zm00028ab.1.gff3.gz    Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3.gz    Zm-P39-REFERENCE-NAM-1.0_Zm00040ab.1.gff3.gz           Zx-mexicana-REFERENCE-YAN-1.0.gff3.gz
Zm-B97-REFERENCE-NAM-1.0_Zm00018ab.1.gff3.gz     Zm-CML52-REFERENCE-NAM-1.0_Zm00019ab.1.gff3.gz   Zm-K0326Y-REFERENCE-SIPPE-1.0_Zm00054a.1.gff3.gz  Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034ab.1.gff3.gz  Zm-PE0075-REFERENCE-TUM-1.0_Zm00017a.1.gff3.gz         Zx-PI566673-REFERENCE-YAN-1.0_Zx00001a.1.gff3.gz
Zm-CML103-REFERENCE-NAM-1.0_Zm00021ab.1.gff3.gz  Zm-CML69-REFERENCE-NAM-1.0_Zm00020ab.1.gff3.gz   Zm-Ki11-REFERENCE-NAM-1.0_Zm00030ab.1.gff3.gz     Zm-Ms71-REFERENCE-NAM-1.0_Zm00035ab.1.gff3.gz   Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0_Zm00008a.1.gff3.gz
Zm-CML228-REFERENCE-NAM-1.0_Zm00022ab.1.gff3.gz  Zm-DK105-REFERENCE-TUM-1.0_Zm00016a.1.gff3.gz    Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1.gff3.gz      Zm-NC350-REFERENCE-NAM-1.0_Zm00036ab.1.gff3.gz  Zm-SK-REFERENCE-YAN-1.0.gff3.gz
Zm-CML247-REFERENCE-NAM-1.0_Zm00023ab.1.gff3.gz  Zm-EP1-REFERENCE-TUM-1.0_Zm00010a.1.gff3.gz      Zm-Ky21-REFERENCE-NAM-1.0_Zm00031ab.1.gff3.gz     Zm-NC358-REFERENCE-NAM-1.0_Zm00037ab.1.gff3.gz  Zm-Tx303-REFERENCE-NAM-1.0_Zm00041ab.1.gff3.gz
Zm-CML277-REFERENCE-NAM-1.0_Zm00024ab.1.gff3.gz  Zm-F7-REFERENCE-TUM-1.0_Zm00011a.1.gff3.gz       Zm-M162W-REFERENCE-NAM-1.0_Zm00033ab.1.gff3.gz    Zm-Oh43-REFERENCE-NAM-1.0_Zm00039ab.1.gff3.gz   Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042ab.1.gff3.gz
Zm-CML322-REFERENCE-NAM-1.0_Zm00025ab.1.gff3.gz  Zm-HP301-REFERENCE-NAM-1.0_Zm00027ab.1.gff3.gz   Zm-M37W-REFERENCE-NAM-1.0_Zm00032ab.1.gff3.gz     Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038ab.1.gff3.gz   Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3.gz
 ```
 - Prepare gene info file
```
$ cat info
Zm00001eb293780 	Zm00001eb.1 	Zm-B73-REFERENCE-NAM-5.0 	MaizeGDB 
ZMex06.022965	.	Zx-mexicana-REFERENCE-YAN-1.0	yan
Zm00015a028623	.	Zm-SK-REFERENCE-YAN-1.0	yan
Zm00004b031071 	Zm00004b.1 	Zm-W22-REFERENCE-NRGENE-2.0 	MaizeGDB 
Zm00010a027938 	Zm00010a.1 	Zm-EP1-REFERENCE-TUM-1.0 	MaizeGDB 
Zm00011a028100 	Zm00011a.1 	Zm-F7-REFERENCE-TUM-1.0 	MaizeGDB 
Zm00014a015747 	Zm00014a.1 	Zm-Mo17-REFERENCE-CAU-1.0 	MaizeGDB 
Zm00018ab307700 	Zm00018ab.1 	Zm-B97-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00019ab281760 	Zm00019ab.1 	Zm-CML52-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00020ab297860 	Zm00020ab.1 	Zm-CML69-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00021ab298350 	Zm00021ab.1 	Zm-CML103-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00022ab298760 	Zm00022ab.1 	Zm-CML228-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00023ab302620 	Zm00023ab.1 	Zm-CML247-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00024ab299240 	Zm00024ab.1 	Zm-CML277-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00025ab304920 	Zm00025ab.1 	Zm-CML322-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00026ab296980 	Zm00026ab.1 	Zm-CML333-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00027ab299040 	Zm00027ab.1 	Zm-HP301-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00028ab301630 	Zm00028ab.1 	Zm-Il14H-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00029ab307720 	Zm00029ab.1 	Zm-Ki3-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00030ab294810 	Zm00030ab.1 	Zm-Ki11-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00031ab305130 	Zm00031ab.1 	Zm-Ky21-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00032ab306600 	Zm00032ab.1 	Zm-M37W-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00033ab311410 	Zm00033ab.1 	Zm-M162W-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00034ab314350 	Zm00034ab.1 	Zm-Mo18W-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00035ab305060 	Zm00035ab.1 	Zm-Ms71-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00036ab302870 	Zm00036ab.1 	Zm-NC350-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00037ab298070 	Zm00037ab.1 	Zm-NC358-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00038ab300540 	Zm00038ab.1 	Zm-Oh7B-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00039ab299570 	Zm00039ab.1 	Zm-Oh43-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00040ab312700 	Zm00040ab.1 	Zm-P39-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00041ab305300 	Zm00041ab.1 	Zm-Tx303-REFERENCE-NAM-1.0 	MaizeGDB 
Zm00042ab307090 	Zm00042ab.1 	Zm-Tzi8-REFERENCE-NAM-1.0 	MaizeGDB
```
You can use the script `fetch_maize_pangene.py` to get this information automatically (input: B73 v5 gene id; output: gene info file):
```
python3 fetch_maize_pangene.py Zm00001eb293780 info
```


