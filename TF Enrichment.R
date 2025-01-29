library(TFBSTools)
library(JASPAR2020)


pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

Control <- AddMotifs(Control, genome = BSgenome.Drerio.UCSC.danRer11, pfm = pfm)

library(chromVAR)
library(presto)
Control <- RunChromVAR(Control, genome = BSgenome.Drerio.UCSC.danRer11)
DefaultAssay(Control) <- "chromvar"
DA_motifs_ct <- wilcoxauc(Control, group_by = "celltype.atac", seurat_assay = "chromvar") %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol),
                           df_pfm$id)[feature])

enriched_motifs_ct <- DA_motifs_ct %>%
  filter(padj < 0.01 & auc > 0.7) %>%
  group_by(group)

top_motifs_ct <- top_n(enriched_motifs_ct, 3, wt=auc)
bluered_colscheme <- colorRampPalette(rev(c("red", "grey", "lightblue")))
FeaturePlot(Control,
            features = "tfap2a",
            cols = bluered_colscheme(30),
            reduction = "umap.atac",
            ncol = 3) 
