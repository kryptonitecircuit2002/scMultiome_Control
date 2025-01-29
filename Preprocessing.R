library(Seurat)
library(Signac)
library(BSgenome.Drerio.UCSC.danRer11)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)

sampleID = "Control_cerbx_script"
setwd("D:/Transit/Single_Cell+ATAC/Control/ATAC/Data+Analysis/Analysis.Control/Cerebrocortex")

#Read files
input.10X.Control <- Read10X_h5("D:/Transit/Single_Cell+ATAC/Combined_Analysis/Control_sample/filtered_feature_bc_matrix_control.h5")
fragpath.Control <- "D:/Transit/Single_Cell+ATAC/Combined_Analysis/Control_sample/atac_fragments_control.tsv.gz"
metadata.Control <- read.csv(
  file = "D:/Transit/Single_Cell+ATAC/Combined_Analysis/Control_sample/per_barcode_metrics_control.csv",
  header = TRUE,
  row.names = 1
)
"D:\Transit\Single_Cell+ATAC\Control\ATAC\Data+Analysis\Analysis.Control\Danio_rerio.GRCz11.105.filtered.gtf"
## Create Seurat object
Control <- CreateSeuratObject(counts = input.10X.Control$`Gene Expression`, 
                           assay = 'RNA',
                           project = sampleID,
                           meta.data = metadata.Control)

Control[["percent.mt"]] <- PercentageFeatureSet(Control, pattern = "^MT-")

#Get annotations 
#annotaion from GTF
gref.path = "D:/Transit/Single_Cell+ATAC/Control/ATAC/Data+Analysis/Analysis.Control/Danio_rerio.GRCz11.105.filtered.gtf"
gene.coords.zf <- rtracklayer::import(gref.path)
genome(gene.coords.zf) <- 'GRCz11'
gene.coords.zf <- gene.coords.zf[! is.na(gene.coords.zf$gene_name),]


##Add ATAC data
atac_counts <- input.10X.Control$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

Control[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'GRCz11',
  fragments = fragpath.Control,
  min.cells = 10,
  annotation = gene.coords.zf
)

#QC
DefaultAssay(Control) <- "ATAC"
Control <- NucleosomeSignal(Control)
Control <- TSSEnrichment(Control, fast = FALSE)
Control$pct_reads_in_peaks <- Control$atac_peak_region_fragments / Control$atac_fragments * 100

#Peak calling
system('bash -c "/home/ayush/.local/bin/macs2 callpeak -t /mnt/d/Transit/Single_Cell+ATAC/Combined_Analysis/Control_sample/atac_fragments_control.tsv.gz -g 1.4e+09 -f BED --nomodel --extsize 200 --shift -100 -n Control_cerbx_script --outdir /mnt/d/Transit/Single_Cell+ATAC/Control/ATAC/Data+Analysis/Analysis.Control/Cerebrocortex" ', 
                wait = TRUE,  ignore.stderr = FALSE,  ignore.stdout = FALSE)
peak_path <- "D:/Transit/Single_Cell+ATAC/Control/ATAC/Data+Analysis/Analysis.Control/Cerebrocortex/Control_cerbx_script_peaks.narrowPeak"
peaks <- rtracklayer::import(peak_path, format = "narrowPeak")
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = 
                          invert = TRUE)

# quantify counts in each peak
macs_count <- FeatureMatrix(fragments = Fragments(Control),
                            features = peaks,
                            cells = colnames(Control))

Control[['ATAC']] <- CreateChromatinAssay(
  counts = macs_count,
  sep = c(":", "-"),
  genome = 'GRCz11',
  fragments = fragpath.Control,
  min.cells = 10,
  annotation = gene.coords.zf
)

DensityScatter(Control, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
VlnPlot(
  object = Control,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

Control <- subset(
  x = Control,
  subset = nCount_ATAC < 50000 &
    nCount_RNA < 10000 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 1 &
    TSS.enrichment < 10)
