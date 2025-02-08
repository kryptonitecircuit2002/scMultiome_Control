#Find DEGs and DARs
DefaultAssay(Control) <- "SCT"
Idents(Control) <- Control$celltype.rna
DEGs <- FindAllMarkers(Control,
                       group.by = "celltype.rna",
                       min.pct = 0.1,
                       min.diff.pct = 0.05,
                       only.pos = T)

DefaultAssay(Control) <- "ATAC"
Idents(Control) <- Control$celltype.atac
da_peaks <- FindAllMarkers(Control, 
                           group.by = "celltype.atac",
                           min.pct = 0.1,
                           min.diff.pct = 0.05,
                           only.pos = T)

#Peak-gene association
library(BSgenome.Drerio.UCSC.danRer11)
Control <- RegionStats(Control, genome = BSgenome.Drerio.UCSC.danRer11)
set.seed(2021)
Control <- LinkPeaks(
  Control, 
  peak.assay = 'ATAC', 
  expression.assay = 'SCT', 
  pvalue_cutoff = 1,
  score_cutoff = 0,
  method = 'spearman'
)
#include adjusted p-values
library(scales)
library(future)
library(data.table)
links <- Control[['ATAC']]@links
pvalue_adjusted = p.adjust(links$pvalue, method = 'fdr')
links$pvalue.fdr = pvalue_adjusted

#Filtering links
#1. using thresholds on correlation coefficient and p-value
links <- links[intersect(which(links$pvalue.fdr < 0.05), which(abs(links$score) > 0.1)),]

#define CollapseToLongest 
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
  ]
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}
#
gene.ranges <- CollapseToLongestTranscript(ranges = Annotation(Control))
peaks <- StringToGRanges(links$peak, sep = c("-", "-"))

######## Quantify DORC scores for genes in links ########
## resulting in a DORC x cell matrix 

## Use 5 peaks per gene as cutoff based on elbow method
peaks.per.gene <- table(links$gene)
genes = names(which(peaks.per.gene >= 5))   

## Aggregate normalized ATAC-seq counts in all significantly associated peaks per gene
atac <- Control[['ATAC']]@data

DORC <- matrix(
  nrow = length(genes), 
  ncol = ncol(atac),
  dimnames = list(genes, colnames(atac))
)

for (i in 1:length(genes)) {
   peaks <- links$peak[which(links$gene == genes[i])]
 
    if (length(peaks) == 0) { # If no peaks are found for the gene, assign a row of zeros (or some default value) in DORC
     DORC[i,] <- rep(0, ncol(atac))
      } else if (length(peaks) > 1) { # If multiple peaks are found, sum their corresponding rows in atac
           DORC[i,] <- Matrix::colSums(atac[peaks, ])
      } else { #If only one peak is found, just copy its row from atac
         DORC[i,] <- atac[peaks, ]
       }
}

######## Plotting statistics of peak-gene links ########

pal <- pal_aaas()(5)

## Distance to TSS
distToTSS = rep(NA, length = length(links))
for (i in 1:length(distToTSS)) {
  
  peak <- peaks[i]
  gene <- gene.ranges[which(gene.ranges$gene_name == links$gene[i])]
  if (length(gene) == 0) next 
  tss <- resize(x = gene, width = 1, fix = 'start')
  
  if(is.na(precede(peak, gene)))
    distToTSS[i] <- -distance(peak, tss)
  else
    distToTSS[i] <- distance(peak,tss)
  
}

## Histogram of distance to TSS (kb)
df = as.data.frame(links@elementMetadata@listData)
df$distance = distToTSS/1000
df <- df[df$distance > -1000 & df$distance_kb < 1000, ]
p1_new <- ggplot(df, aes(x = distance)) + 
      geom_histogram(bins = 30, color = 'black', fill = pal[5]) +
     xlab('Distance to TSS (kb)') +
     ylab('Count') +
     theme_classic() +
     theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
           axis.title = element_text(face = "bold")
        ) +
    ggtitle("Distribution of Distances to TSS")
p1_new

## Identify genes 'skipped' by peak in each link
fo <- findOverlaps(links, gene.ranges)
fo <- as.data.frame(fo)
fo$gene <- gene.ranges$gene_name[fo$subjectHits]
skipped.genes <- vector('list', length = length(links))

for (i in 1:length(links)) {
  
  idx <- which(fo$queryHits == i)
  g <- setdiff(fo$gene[idx], links$gene[i])
  skipped.genes[[i]] <- g
  
}

skipped.genes.num <- sapply(skipped.genes, function(x) length(unique(x)))

## Histograms of #peaks per gene, #genes per peaks, #genes skipped by peak
peaks.per.gene <- data.frame(ppg = table(links$gene))
genes.per.peak <- data.frame(gpp = table(links$peak))
genes.skipped <- data.frame(skip = skipped.genes.num)
mean(peaks.per.gene$ppg.Freq)
quantile(peaks.per.gene$ppg.Freq, 0.95)
p2_1 = ggplot(peaks.per.gene, aes(x = ppg.Freq)) + 
  geom_histogram(binwidth = 1, color = 'black', fill = pal[1]) + 
  xlab('# Peaks mapped to genes') +
  annotate(geom = 'text', label = 'Mean=2.407\n95%=7', x = 20, y = 1500) +
  theme_classic()

mean(genes.per.peak$gpp.Freq)
sum(genes.per.peak$gpp.Freq == 1) / length(genes.per.peak$gpp.Freq) * 100
p2_2 = ggplot(genes.per.peak, aes(x = gpp.Freq)) + 
  geom_histogram(binwidth = 1, color = 'black', fill = pal[2]) + 
  xlab('# Genes mapped to peaks') +
  scale_x_continuous(breaks = seq(0,10,2)) +
  annotate(geom = 'text', label = 'Mean=1.23\nTo 1 Gene=86.9%', x = 7, y = 4500) +
  theme_classic() +
  theme(axis.title.y = element_blank())

mean(genes.skipped$skip)
sum(genes.skipped$skip == 0) / length(genes.skipped$skip) * 100
p2_3 = ggplot(genes.skipped, aes(x = skip)) + 
  geom_histogram(binwidth = 2, color = 'black', fill = pal[3]) + 
  xlab('# Genes skipped by peak') +
  theme(axis.title.y = element_blank()) +
  annotate(geom = 'text', label = 'Mean=6.07\nNearest Gene=8.51%', x = 25, y = 3500) + 
  theme_classic() +
  theme(axis.title.y = element_blank())

p2_1 + p2_2 + p2_3
p2 = gridExtra::grid.arrange(p2_1, p2_2, p2_3, ncol = 3)

## Plot number of correlated peaks in increasing order vs. gene ranks
peaks.per.gene <- table(links$gene)

y = as.numeric(sort(peaks.per.gene))
df_p3 <- data.frame(ind = 1:length(y), count = y, gene = names(sort(peaks.per.gene)))
p3 = ggplot(df_p3, aes(x = ind, y = count)) +
  geom_line(color = pal[1])+
  geom_point(color = pal[1], size = 0.6) +
  geom_hline(yintercept = 5, linetype = 'longdash', color = 'grey') +
  geom_vline(xintercept = 2718, linetype = 'longdash', color = 'grey') +
  ggrepel::geom_text_repel(aes(label = ifelse(count > 10, gene, '')), size = 2.5, max.overlaps = 16, direction = 'both') +
  theme_classic() +
  xlab('Ranked gene list') +
  ylab('Number of correlated peaks') +
  xlim(0,7000) +
  scale_y_continuous(breaks=c(5,10, 15, 20, 25))
p3

## Plot positive vs. negative correlations
p4 = ggplot(links,aes(x = score)) + 
  geom_histogram(position = 'identity', fill = 'steelblue', color = 'black') + 
  theme_classic() + 
  xlab('Spearman correlation coefficient') +
  ylab('# peak-gene links')
p4
