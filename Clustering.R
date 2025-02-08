#Filter peaks/genes 
set.seed(1234)  
tmp.rna <- Matrix::rowSums(Control[["RNA"]]@layers$counts > 0)
Control[["RNA"]] <- subset(Control[["RNA"]], features = names(which(tmp.rna >= 10)))
tmp.atac <- Matrix::rowSums(Control[["ATAC"]]@counts > 0)
Control[["ATAC"]] <- subset(Control[["ATAC"]], features = names(which(tmp.atac >= 10)))
  
## Normalization, dimensional reduction, and clustering on RNA-seq and ATAC-seq separately
  
# RNA-seq
DefaultAssay(Control) <- 'RNA'
Control <- SCTransform(Control) %>%
RunPCA(npcs = 30, verbose = F) %>%
FindNeighbors(reduction = 'pca', dims = 1:30) %>%
FindClusters(graph.name = 'SCT_snn', algorithm = 3, resolution = 0.5, verbose = FALSE) 
Control <- RunUMAP(Control, reduction = 'pca', dims = 1:30,
                     reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',spread = 0.28)
DimPlot(Control, reduction = "umap.rna", label = T)  

# ATAC-seq
DefaultAssay(Control) <- "ATAC"
Control <- RunTFIDF(Control, method = 3)
Control <- FindTopFeatures(Control, min.cutoff = 'q75')
Control <- RunSVD(Control)
Control <- FindNeighbors(Control, reduction = 'lsi', dims = 2:10, assay = 'ATAC')
Control <- FindClusters(Control, graph.name = 'ATAC_snn', algorithm = 2, resolution = 0.5)
Control <- RunUMAP(Control, reduction = 'lsi', dims = 2:10, assay = 'ATAC',
                    reduction.name = "umap.atac", reduction.key = "atacUMAP_", spread = 0.28)
DimPlot(Control, reduction = "umap.atac", label = T)

# Weighted nearest neighbor (WNN) analysis using both modalities
Control <- FindMultiModalNeighbors(Control,
                                    reduction.list = list("pca", "lsi"),
                                    dims.list = list(1:30, 2:10),
                                    modality.weight.name = 'RNA.weight')
Control <- FindClusters(Control, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
Control <- RunUMAP(Control, nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
DimPlot(Control, reduction = "wnn.umap", label = T)

#Naming clusters
## WNN-derived
celltype.wnn <- rep(NA, length = ncol(Control))
Idents(Control) <- Control$wsnn_res.0.5

celltype.wnn[which(Idents(Control) %in% c(0,10))] <- 'twist 1a+ progenitor'
celltype.wnn[which(Idents(Control) %in% c(1,4))] <- 'mesenchymal'
celltype.wnn[which(Idents(Control) %in% c(3))] <- 'glial/neural'
celltype.wnn[which(Idents(Control) %in% c(2))] <- 'MIX+'
celltype.wnn[which(Idents(Control) %in% c(7))] <- 'foxd3 progenitor 2'
celltype.wnn[which(Idents(Control) %in% c(5))] <- 'neuronal'
celltype.wnn[which(Idents(Control) %in% c(6))] <- 'foxd3 progenitor 1'
celltype.wnn[which(Idents(Control) %in% c(8))] <- 'otic'
celltype.wnn[which(Idents(Control) %in% c(9))] <- 'melanoblast'
celltype.wnn[which(Idents(Control) %in% c(11))] <- 'iridoblast'
celltype.wnn <- factor(celltype.wnn, 
                   levels = c('twist 1a+ progenitor', 'mesenchymal', 'foxd3 progenitor 1', 'foxd3 progenitor 2', 'MIX+', 
                              'melanoblast', 'iridoblast', 'glial/neural',
                              'neuronal', 'otic'), 
                   ordered = T)
Control$celltype.wnn <- celltype.wnn
#RNA Derived
celltype.rna <- rep(NA, length = ncol(Control))
Idents(Control) <- Control$SCT_snn_res.0.5

celltype.rna[which(Idents(Control) %in% c(0,10))] <- 'twist 1a+ progenitor'
celltype.rna[which(Idents(Control) %in% c(1,4))] <- 'mesenchymal'
celltype.rna[which(Idents(Control) %in% c(2))] <- 'glial/neural'
celltype.rna[which(Idents(Control) %in% c(3))] <- 'MIX+'
celltype.rna[which(Idents(Control) %in% c(7))] <- 'foxd3 progenitor 1'
celltype.rna[which(Idents(Control) %in% c(5))] <- 'foxd3 progenitor 2'
celltype.rna[which(Idents(Control) %in% c(6,13))] <- 'neuronal'
celltype.rna[which(Idents(Control) %in% c(8))] <- 'otic'
celltype.rna[which(Idents(Control) %in% c(9))] <- 'melanoblast'
celltype.rna[which(Idents(Control) %in% c(11))] <- 'iridoblast'
celltype.rna[which(Idents(Control) %in% c(12))] <- 'skeletal muscle progenitor'
celltype.rna <- factor(celltype.rna, 
                   levels = c('twist 1a+ progenitor', 'mesenchymal', 'foxd3 progenitor 1', 'foxd3 progenitor 2', 'MIX+', 
                              'melanoblast', 'iridoblast', 'glial/neural',
                              'neuronal', 'otic', 'skeletal muscle progenitor'), 
                   ordered = T)

Control$celltype.rna <- celltype.rna
#ATAC derived
gene.activities_control <- GeneActivity(Control)

celltype.atac <- rep(NA, length = ncol(Control))
Idents(Control) <- Control$ATAC_snn_res.0.5

celltype.atac[which(Idents(Control) %in% c(0,9))] <- 'twist 1a+ progenitor'
celltype.atac[which(Idents(Control) %in% c(1,6))] <- 'mesenchymal'
celltype.atac[which(Idents(Control) %in% c(4))] <- 'glial/neural'
celltype.atac[which(Idents(Control) %in% c(5))] <- 'MIX+'
celltype.atac[which(Idents(Control) %in% c(3))] <- 'foxd3 progenitor'
celltype.atac[which(Idents(Control) %in% c(7))] <- 'neuronal'
celltype.atac[which(Idents(Control) %in% c(8))] <- 'otic'
celltype.atac[which(Idents(Control) %in% c(10))] <- 'melanoblast'
celltype.atac[which(Idents(Control) %in% c(2))] <- 'iridoblast'

celltype.atac <- factor(celltype.atac, 
                       levels = c('twist 1a+ progenitor', 'mesenchymal', 'foxd3 progenitor', 'MIX+', 
                                  'melanoblast', 'iridoblast', 'glial/neural',
                                  'neuronal', 'otic'), 
                       ordered = T)
Control$celltype.atac <- celltype.atac

#compute RNA and ATAC cluster concordance
contingency_table <- table(Control$celltype.rna, Control$celltype.atac)
melted_table <- melt(contingency_table)
ggplot(melted_table, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  scale_fill_gradientn(colors = c("white", "orange", "red")) +
  labs(title = "RNA-ATAC Concordance", 
       x = "RNA annotations", 
       y = "ATAC annotations", 
       fill = "Count") + 
  theme_bw()

#compute geneactivity
DefaultAssay(Control) <- "ATAC"
gene.activities <- GeneActivity(Control, assay = "ATAC")
Control[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
Control <- NormalizeData(Control, assay = 'GeneActivity')

DefaultAssay(Control) <- "GeneActivity"
Idents(Control) <- Control$celltype.atac
g1 <- DotPlot(Control, features = c("twist1a", "grem2b", "col11a1b", "lamc3", "foxd3", "zeb2a", "sox10", "crestin", "aox5", "tfec", "paics", "mitfa", "dct", "kita", "ltk", "alx4a", "gfap", "pou3f1", "elavl3", "elavl4", "neurod1", "oc90", "epcam", "tpma", "musk")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() +
  coord_flip() +
  ggtitle("Gene Activity") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

DefaultAssay(Control) <- "SCT"
Idents(Control) <- Control$celltype.rna
r1 <- DotPlot(Control, features = c("twist1a", "grem2b", "col11a1b", "lamc3", "foxd3", "zeb2a", "sox10", "crestin", "aox5", "tfec", "paics", "mitfa", "dct", "kita", "ltk", "alx4a", "gfap", "pou3f1", "elavl3", "elavl4", "neurod1", "oc90", "epcam", "tpma", "musk")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() +
  coord_flip() +
  ggtitle("Gene Expression") +
  theme(axis.title.y = element_blank())

r1 + g1 + plot_layout(guides = "collect") & theme(legend.position = "right")
