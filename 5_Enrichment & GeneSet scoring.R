rm(list = ls())
options(stringsAsFactors = F)

library(Seurat)
library(tidyverse)
library(ggsci)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GSEABase)
library(RColorBrewer)

library(zellkonverter)
sce <- readH5AD("C:/path/to/apc_forSeurat.h5ad", verbose = TRUE)
apc <- as.Seurat(sce, counts = "X", data = NULL)
save(apc, file='apc_from_scanpy.Rdata')
# load('apc_from_scanpy.Rdata')

apc <- FindVariableFeatures(apc, selection.method = "vst", nfeatures = 3000)
apc <- ScaleData(apc)
Idents(apc) <- 'cluster'

#### GPCR list ####
gpcr_human <- read.csv('input/GPCR-new_human.csv', header = F)
library(homologene)
gpcr_mouse <- human2mouse(gpcr_human$V1)
gpcr_mouse <- c(gpcr_mouse$mouseGene, str_to_title(gpcr_human$V1))
gpcr_mouse <- gpcr_mouse[!duplicated(gpcr_mouse)]
gpcr_mouse <- gpcr_mouse[gpcr_mouse %in% rownames(apc)]
write.csv(gpcr_mouse, file='gpcr_Nat_mouse.csv')


#### GO:BP enrichment: driver genes ####==============================================
driver <- read.csv(paste0(scanpyPath, 'drivers.csv'), header = T)
driver <- driver[,1:3]
colnames(driver) <- c('symbol', 'corr', 'pval')
driver <- driver %>% filter(pval<=0.05)

drivergene <- driver$symbol[1:500]
drivergene_r <- rev(driver$symbol)[1:500]

driver_ego <- data.frame(enrichGO(gene=c(drivergene,drivergene_r), OrgDb='org.Mm.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff=1, qvalueCutoff=1))
driver_ego <- driver_ego %>% 
  dplyr::filter(p.adjust<0.05) %>% 
  dplyr::mutate(log10Padj=-log10(p.adjust), regulated='up') %>% 
  dplyr::select(Description, p.adjust, log10Padj, regulated)

driver_path <- c(
  'extracellular matrix organization',
  'cell-substrate adhesion',
  'small GTPase mediated signal transduction',
  'ERK1 and ERK2 cascade',
  'negative regulation of cell motility',
  'regulation of angiogenesis',
  'connective tissue development',
  'Ras protein signal transduction',
  'regulation of cell morphogenesis',
  'cellular response to peptide',
  'positive regulation of lipid metabolic process',
  'Wnt signaling pathway',
  'developmental cell growth',
  'fat cell differentiation',
  'cellular response to BMP stimulus'
)

driver_ego <- driver_ego[driver_ego$Description %in% driver_path,]

ggplot(driver_ego)+
  geom_col(aes(x = reorder(Description, log10Padj), y = log10Padj),
           fill = '#438CBF', width = 0.6) +
  geom_hline(yintercept = -log10(0.05), color = 'Gray50', linetype='longdash')+
  coord_flip()+
  scale_y_continuous(name = "-log10(p adjust)", 
                     limits = c(0,45), breaks = c(0,10,20,30,40), expand = c(0,0)) +
  theme_bw()+
  theme(text = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank()
  )
ggsave('GOBPenrich_driver.pdf', units = 'in', width = 9, height = 8.5)



#### AUCell ####==================================================================
library(AUCell)
library(GSEABase)

# gene expression matrix
geneMatrix4aucell <- apc[["originalexp"]]@data
adipogenic_score <- c('Pparg','Cebpa','Cebpb','Bmp2','Bmp4','Smad4','Smad1','Smad5','Smad8','Zfp423',
                      'Zfp467','Zfp217','Klf4','Klf5','Klf13','Klf15','Egr2','Creb5','Bcl6','Ebf1',
                      'Rxra','Nr3c1','Nfil3','Igf1','Ffar4')
geneSets4aucell <- c(
  GeneSet(adipogenic_score, setName = 'adipogenic_score')
)
geneSets4aucell <- GeneSetCollection(geneSets4aucell)
geneSets4aucell <- subsetGeneSets(geneSets4aucell, rownames(apc))
names(geneSets4aucell)

# run AUCell
cells_AUC <- AUCell_run(geneMatrix4aucell, geneSets4aucell, 
                        keepZeroesAsNA = TRUE,
                        aucMaxRank = ceiling(0.08 * nrow(geneMatrix4aucell)))

# add results to metadata
View(apc@meta.data)
aucMatrix <- t(getAUC(cells_AUC))
aucMatrix <- as.data.frame(aucMatrix)
identical(rownames(apc@meta.data), rownames(aucMatrix))
apc@meta.data <- cbind(apc@meta.data, aucMatrix)


#### GPCR expression ####=============================================
gpcr <- read.csv('output/gpcr_Nat_mouse.csv')
gpcr <- gpcr$symbol_mouse
gpcr_mean <- AverageExpression(apc, features = gpcr)
gpcr_mean <- as.data.frame(gpcr_mean$originalexp)
gpcr_mean <- gpcr_mean[,1:2]
gpcr_mean$mean <- rowMeans(gpcr_mean)
gpcr_mean <- gpcr_mean %>% arrange(desc(mean))
gpcr_mean <- gpcr_mean[rowSums(gpcr_mean)!=0,]
write.csv(gpcr_mean, file = 'output/gpcr_expression_2clusters.csv')
gpcr_mean <- gpcr_mean %>% rownames_to_column(var='symbol')
master_gpcr <- read.csv('input/master_driver-gpcr.csv')
gpcr_mean$master = ifelse(gpcr_mean$symbol %in% master_gpcr$symbol, 'master', 'ordinary')

ggplot(gpcr_mean)+
  geom_point(aes(x=reorder(symbol, mean), y=mean, color=master, alpha=master), size=2.5, shape=16)+
  scale_color_manual(values = c('#d62728','gray50'))+
  scale_alpha_manual(values = c(0.9,0.5))+
  scale_y_continuous(name = "Normalized expression", 
                     limits = c(-1,13.5), breaks = c(0,2,4,6,8,10,12),expand = c(0,0)) +
  theme_classic() + theme(panel.grid = element_blank())+
  theme(text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14))
ggsave('gpcr_exp_point.pdf', units = 'in', width = 8.5, height = 5)

