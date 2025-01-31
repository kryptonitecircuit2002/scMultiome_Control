# scMultiome_Control

Samples: 24 hpf embryos from zebrafish sox10+:GFP transgenic line \
Kit used: 10X 3' GEX + ATAC scMultiome \
Genome: GRCz11 ver 105 (custom made using cellranger mkref)

Files Description:\
```Preprocessing.R```: create seurat objects and perform QC \
```Clustering.R```: Perform RNA, ATAC and Wnn clustering \
```Peak-gene association.R```: Link peaks and genes and perform statistical analysis \
```TF enrichment.R```: Transcription Factor Enrichment analysis
