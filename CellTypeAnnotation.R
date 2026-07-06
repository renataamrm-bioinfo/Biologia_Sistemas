# ==============================================================================
# PIPELINE PRINCIPAL: LOOP DE QC DINÂMICO E MARCADORES GENÔMICOS
# ==============================================================================

# ------------------------------------------------------------------------------
# ETAPA 0: Carregar Pacotes e Configurar Ambiente
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scuttle)
library(scDblFinder)
library(ggplot2)
library(SeuratDisk)
library(sctransform)
library(writexl)

# Garantir a reprodutibilidade computacional (mesmos resultados em algoritmos estocásticos)
set.seed(1234)

# Aumentar o limite de memória para lidar com múltiplos datasets massivos (alocação de 15 GB)
# Se houver erros de memória, considere aumentar este valor se sua máquina permitir.
options(future.globals.maxSize = 15 * 1024^3)

# ==============================================================================
# ORGANIZAÇÃO DE DIRETÓRIOS (INPUTS E OUTPUTS)
# ==============================================================================
# Defina o diretório base aqui. Todo o resto do script se adaptará automaticamente.
base_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/"

# Subdiretórios padronizados para organizar o fluxo de trabalho
data_dir   <- base_dir                       # Onde estão as matrizes raw (PT1, PT2, etc.)
plot_dir   <- paste0(base_dir, "plots/")     # Destino para todos os gráficos salvos
result_dir <- paste0(base_dir, "results/")   # Destino para tabelas excel e objetos RDS

# Criar os diretórios de saída caso eles não existam no seu computador
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# ETAPA 1: Carregar Matrizes scRNA-seq (10X Genomics)
# ------------------------------------------------------------------------------
print("Carregando os diretórios de matrizes do 10X Genomics para a memória...")
PT1 <- Read10X(data.dir = paste0(data_dir, "PT1"))
PT2 <- Read10X(data.dir = paste0(data_dir, "PT2"))
PT3 <- Read10X(data.dir = paste0(data_dir, "PT3"))
RT1 <- Read10X(data.dir = paste0(data_dir, "RT1"))
RT2 <- Read10X(data.dir = paste0(data_dir, "RT2"))
RT3 <- Read10X(data.dir = paste0(data_dir, "RT3"))
NT1 <- Read10X(data.dir = paste0(data_dir, "NT1"))
NT2 <- Read10X(data.dir = paste0(data_dir, "NT2"))

# ------------------------------------------------------------------------------
# ETAPA 2: Inicializar Objetos Seurat Iniciais
# ------------------------------------------------------------------------------
# Filtramos genes expressos em menos de 3 células e células com menos de 200 genes
PT1.seu <- CreateSeuratObject(counts = PT1, project = "PT1", min.cells = 3, min.features = 200)
PT2.seu <- CreateSeuratObject(counts = PT2, project = "PT2", min.cells = 3, min.features = 200)
PT3.seu <- CreateSeuratObject(counts = PT3, project = "PT3", min.cells = 3, min.features = 200)
RT1.seu <- CreateSeuratObject(counts = RT1, project = "RT1", min.cells = 3, min.features = 200)
RT2.seu <- CreateSeuratObject(counts = RT2, project = "RT2", min.cells = 3, min.features = 200)
RT3.seu <- CreateSeuratObject(counts = RT3, project = "RT3", min.cells = 3, min.features = 200)
NT1.seu <- CreateSeuratObject(counts = NT1, project = "NT1", min.cells = 3, min.features = 200)
NT2.seu <- CreateSeuratObject(counts = NT2, project = "NT2", min.cells = 3, min.features = 200)

# Criar uma lista para iterar sobre todas as amostras
seurat_list <- c('PT1' = PT1.seu, 'PT2' = PT2.seu, 'PT3' = PT3.seu,
                 'RT1' = RT1.seu, 'RT2' = RT2.seu, 'RT3' = RT3.seu,
                 'NT1' = NT1.seu, 'NT2' = NT2.seu)

# ------------------------------------------------------------------------------
# ETAPA 3: Funções Auxiliares para Controle de Qualidade (QC)
# ------------------------------------------------------------------------------

# Função para encontrar "doublets" (gotículas que capturaram 2 células por engano)
finddoublet <- function(seurat_obj){
  set.seed(34) 
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder(sce_obj)
  # Extraímos a classificação e o escore e devolvemos ao objeto Seurat original
  seurat_obj$scDblFinder.class <- sce_obj$scDblFinder.class
  seurat_obj$scDblFinder.score <- sce_obj$scDblFinder.score
  return(seurat_obj)
}

# Função para calcular Limites Dinâmicos (MAD - Median Absolute Deviation)
# Evita o uso de limiares arbitrários (ex: cortar sempre em 5000 genes)
nMAD <- function(x, nmads=3){
  xm <- median(x)
  md <- median(abs(x-xm))
  mads <- xm + nmads*md
  return(mads)
}

# ------------------------------------------------------------------------------
# ETAPA 4: Loop de Processamento e Filtro de Qualidade (QC)
# ------------------------------------------------------------------------------
nmad = 5 # Quão rigoroso será o limite dinâmico
count = 1

# Vetores para guardar os limiares gerados para cada amostra (útil para os gráficos)
sample.nfeaure.cut <- c()
sample.ncount.cut <- c()
sample.mt.cut <- c()
seurat_list_qc <- c()

for (obj in seurat_list){
  print(names(seurat_list)[count])
  
  # 1. Calcular a porcentagem de genes mitocondriais (células mortas/estressadas)
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # 2. Visualização Pré-QC
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  
  # 3. Detectar doublets
  obj = finddoublet(obj)
  
  # 4. Calcular cortes máximos usando MAD para cada métrica
  nfeature.upcut <- ceiling(nMAD(obj$nFeature_RNA, nmad))
  ncount.upcut <- ceiling(nMAD(obj$nCount_RNA, nmad))
  permt.upcut <- ceiling(nMAD(obj$percent.mt, nmad))
  
  # Guardar esses cortes nos vetores
  sample.nfeaure.cut <- c(sample.nfeaure.cut, nfeature.upcut)
  sample.ncount.cut <- c(sample.ncount.cut, ncount.upcut)
  sample.mt.cut <- c(sample.mt.cut, permt.upcut)
  
  # 5. Aplicar os filtros de fato no objeto Seurat
  # Removemos valores muito altos (outliers) e muito baixos (gotículas vazias/ruído)
  obj.filt <- subset(obj, subset = nFeature_RNA <= nfeature.upcut & nFeature_RNA > nfeature.upcut/20)
  obj.filt <- subset(obj.filt, subset = nCount_RNA <= ncount.upcut & nCount_RNA > ncount.upcut/20)
  
  # Filtro mitocondrial: usamos o limite do MAD ou no máximo 25% (o que for menor)
  obj.filt <- subset(obj.filt, subset = percent.mt <= min(permt.upcut, 25))
  
  # Mantemos apenas as células singulares (singlets)
  obj.filt <- subset(obj.filt, subset = scDblFinder.class %in% c('singlet'))
  
  Idents(obj.filt) <- names(seurat_list)[count]
  
  # Visualização Pós-QC (Sanity check)
  print(VlnPlot(obj.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  
  # Salvar o objeto filtrado na lista final
  seurat_list_qc <- c(seurat_list_qc, obj.filt)
  count = count + 1
}

# Organizar os nomes das listas finais
names(sample.mt.cut) <- names(seurat_list)
names(sample.ncount.cut) <- names(seurat_list)
names(sample.nfeaure.cut) <- names(seurat_list)
names(seurat_list_qc) <- names(seurat_list)


# ------------------------------------------------------------------------------
# ETAPA 5: Preparar Dados para os Gráficos Suplementares (QC)
# ------------------------------------------------------------------------------
samples <- names(seurat_list)

# Dataframes estruturados para desenhar as linhas tracejadas nos gráficos
samples.ncount.cut.summary <- data.frame(
  orig.ident = c(samples, samples), 
  types = c(rep('up.cut', length(samples)), rep('down.cut', length(samples))),
  values = c(sample.ncount.cut, sample.ncount.cut/20)
)

samples.nfeature.cut.summary <- data.frame(
  orig.ident = c(samples, samples), 
  types = c(rep('up.cut', length(samples)), rep('down.cut', length(samples))),
  values = c(sample.nfeaure.cut, sample.nfeaure.cut/20)
)

# Genes mitocondriais não precisam de corte inferior
samples.mt.cut.summary <- data.frame(
  orig.ident = c(samples), 
  types = c(rep('up.cut', length(samples))),
  values = c(sample.mt.cut)
)

# Mesclar os dados brutos (pre-QC) apenas para geração de gráficos comparativos
merge.seu <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
merge.seu[['percent.mt']] <- PercentageFeatureSet(merge.seu, pattern = "^MT-")
metadata <- merge.seu@meta.data

# ==============================================================================
# GERAR E SALVAR GRÁFICOS PRÉ-QC
# ==============================================================================
print("#### Gerando gráficos Pré-QC ####" )

# 1. Gráfico de Barras: Número de células
metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("NCells")
ggsave(paste0(plot_dir, "01_PreQC_NCells_barplot.png"), width = 6, height = 4, dpi = 300)

# 2. Gráfico de Densidade: nCount_RNA
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)
ggsave(paste0(plot_dir, "02_PreQC_nCount_density.png"), width = 10, height = 6, dpi = 300)

# 3. Gráfico de Densidade: nFeature_RNA
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_classic() +
  geom_vline(data = samples.nfeature.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)
ggsave(paste0(plot_dir, "03_PreQC_nFeature_density.png"), width = 10, height = 6, dpi = 300)

# 4. Gráfico de Densidade: percent.mt
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_classic() +
  geom_vline(data = samples.mt.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)
ggsave(paste0(plot_dir, "04_PreQC_percentMT_density.png"), width = 10, height = 6, dpi = 300)

# 5. Gráfico de Dispersão (Scatter): nCount vs nFeature
metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + scale_colour_gradient(low = "gray90", high = "black") +
  geom_smooth(se=TRUE,level=0.9) + scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F)+
  geom_hline(data = samples.nfeature.cut.summary, aes(yintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)
ggsave(paste0(plot_dir, "05_PreQC_Scatter_Count_vs_Feature.png"), width = 12, height = 8, dpi = 300)

# ==============================================================================
# GERAR E SALVAR GRÁFICOS PÓS-QC
# ==============================================================================
print("#### Gerando gráficos Pós-QC ####" )

# Mesclar a lista limpa e extrair metadados para avaliar o sucesso do filtro
merge.qc.seu <- merge(x=seurat_list_qc[[1]], y=seurat_list_qc[2:length(seurat_list_qc)])
metadata.qc <- merge.qc.seu@meta.data

# 1. Gráfico de Barras: Número de células (Pós-QC)
metadata.qc %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave(paste0(plot_dir, "01_PostQC_NCells_barplot.png"), width = 6, height = 4, dpi = 300)

# 2. Gráfico de Densidade: nCount_RNA (Pós-QC)
metadata.qc %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)
ggsave(paste0(plot_dir, "02_PostQC_nCount_density.png"), width = 10, height = 6, dpi = 300)

# 3. Gráfico de Densidade: nFeature_RNA (Pós-QC)
metadata.qc %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_classic() +
  geom_vline(data = samples.nfeature.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)
ggsave(paste0(plot_dir, "03_PostQC_nFeature_density.png"), width = 10, height = 6, dpi = 300)

# 4. Gráfico de Densidade: percent.mt (Pós-QC)
metadata.qc %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_classic() +
  geom_vline(data = samples.mt.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)
ggsave(paste0(plot_dir, "04_PostQC_percentMT_density.png"), width = 10, height = 6, dpi = 300)

# 5. Gráfico de Dispersão: nCount vs nFeature (Pós-QC)
metadata.qc %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + scale_colour_gradient(low = "gray90", high = "black") +
  geom_smooth(se=TRUE,level=0.9) + scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, aes(xintercept = values), linetype = "dashed", show.legend = F)+
  geom_hline(data = samples.nfeature.cut.summary, aes(yintercept = values), linetype = "dashed", show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)
ggsave(paste0(plot_dir, "05_PostQC_Scatter_Count_vs_Feature.png"), width = 12, height = 8, dpi = 300)


# ==============================================================================
# ANÁLISES PRINCIPAIS POR TECIDO (NTs, TTs, PTs, RTs)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. TECIDOS NORMAIS (NTs)
# ------------------------------------------------------------------------------
# Mesclar amostras normais
NTs.seu <- merge(seurat_list_qc$NT1, y = c(seurat_list_qc$NT2), 
                 add.cell.ids = c("NT1", "NT2"), project = "scRNA_NTs")

# Normalização SCTransform (método glmGamPoi regressando genes mitocondriais)
NTs.seu <- SCTransform(NTs.seu, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# Redução de dimensionalidade linear (PCA) e não linear (UMAP/t-SNE)
NTs.seu <- RunPCA(NTs.seu, npcs = 30, verbose = FALSE)
NTs.seu <- RunUMAP(NTs.seu, reduction = "pca", dims = 1:30)
NTs.seu <- RunTSNE(NTs.seu, reduction = "pca", dims = 1:30)

# Construir grafo KNN e definir agrupamentos (clusters)
NTs.seu <- FindNeighbors(NTs.seu, reduction = "pca", dims = 1:30)
NTs.seu <- FindClusters(NTs.seu, resolution = 0.9, algorithm = 2)

DimPlot(NTs.seu, label=TRUE, repel = T)

# Encontrar genes marcadores para cada cluster em comparação aos demais
NTs.seu <- PrepSCTFindMarkers(NTs.seu)
NTs_markers <- FindAllMarkers(object = NTs.seu, assay='SCT', only.pos = TRUE, logfc.threshold = 0.25)

# Selecionar o Top 20 marcadores por cluster (baseado no log2FoldChange)
NTs_top <- NTs_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Salvar tabelas Excel no diretório de resultados
sheetlist <- list("allPosMarkers"=NTs_markers,"Top20Markers"=NTs_top)
write_xlsx(sheetlist, paste0(result_dir, 'NTs_markers.xlsx'))

# Validação Biológica - Marcadores Canônicos
NTs.marker <- c("KRT15","KRT16","KRT17",         # Basais
                "ESRP1","ELF3","RARRES1",        # Luminais
                "TOP2A", "CDK1","MKI67","CENPF", # Progenitores Luminais
                "COL4A6","COL4A5",               # Fibroblastos
                "S100A10", "ID1",                # Células Endoteliais
                "CXCL1", "CXCL8",                # Células Mielóides
                "CD24")                          # Luminal/Epitelial Geral

# Plot de Violino (Stacked) para inspecionar onde os marcadores caíram
VlnPlot(NTs.seu, features = NTs.marker, stack=T, pt.size=0, flip = T, add.noise = T) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), axis.text.x = element_text(colour = 'black',size = 14),
        legend.position = 'none')
ggsave(paste0(plot_dir, "NTs_Canonical_Markers_Violin.png"), width = 10, height = 8, dpi = 300)

# Renomear os clusters com base no perfil genético encontrado
NTs.seu <- RenameIdents(object = NTs.seu,
                        "0" = "Luminal cells", "1" = "Basal cells", "2" = "Endothelial cells", 
                        "3" = "Luminal cells", "4" = "Fibroblasts", "5" = "Basal cells",
                        "6" = "Luminal Progenitor", "7" = "Luminal Progenitor",
                        "8" = "Myeloid cells", "9" = "Luminal Progenitor", "10" = "Luminal cells")

NTs.seu$cell.annot <- Idents(NTs.seu)
DimPlot(NTs.seu, label = FALSE, repel = TRUE, label.size = 6) + ggtitle("Annotated Normal Tissue (Clean Pipeline)")
ggsave(paste0(plot_dir, "Normal Tissue Anoot.png"), width = 12, height = 8, dpi = 300)

# Salvar objeto Seurat no diretório de resultados
saveRDS(NTs.seu, file = paste0(result_dir, "NTs_classified.rds"))


# ------------------------------------------------------------------------------
# 2. TUMORES TOTAIS (TTs) E PRIMÁRIOS (PTs)
# ------------------------------------------------------------------------------
# Preparação dos dados tumorais em blocos
seu.TTs.part1 <- merge(seurat_list_qc$PT1, y = c(seurat_list_qc$PT2, seurat_list_qc$RT1, seurat_list_qc$RT2), 
                       add.cell.ids = c("PT1", "PT2", "RT1", "RT2"), project = "scRNA_09302020")
seu.TTs.part2 <- merge(seurat_list_qc$PT3, y = c(seurat_list_qc$RT3), 
                       add.cell.ids = c("PT3", "RT3"), project = "scRNA_01112021")

# Integração Total 
TTs.combined <- merge(seu.TTs.part1, y = seu.TTs.part2, add.cell.ids = c("Batch1", "Batch2"), project = "TTs_Integrated")

# Correção Seurat V5: Utiliza-se variable.features.n
TTs.combined <- SCTransform(TTs.combined, vst.flavor = "v2", variable.features.n = 3000, verbose = FALSE)
TTs.combined.sct <- TTs.combined
TTs.combined.sct <- RunPCA(TTs.combined.sct, npcs = 30, verbose = FALSE)
TTs.combined.sct <- RunTSNE(TTs.combined.sct, reduction = "pca", dims = 1:30)
TTs.combined.sct <- RunUMAP(TTs.combined.sct, reduction = "pca", dims = 1:30)
TTs.combined.sct <- FindNeighbors(TTs.combined.sct, reduction = "pca", dims = 1:30)
TTs.combined.sct <- FindClusters(TTs.combined.sct, resolution = 0.9)

saveRDS(TTs.combined.sct, file = paste0(result_dir, "TTS.rds"))

DimPlot(TTs.combined.sct, label=TRUE, repel = T)
ggsave(paste0(plot_dir, "TT Tissue Clsuter.png"), width = 12, height = 8, dpi = 300)

# Isolando os Tumores Primários (PTs)
PTs.combined.sct <- subset(TTs.combined.sct, subset = orig.ident %in% c('PT1','PT2','PT3'))
PTs.combined.sct <- RunPCA(PTs.combined.sct, npcs = 30, verbose = FALSE)

# Integração das camadas (batch effect removal) limpa e robusta via CCA
PTs.combined.sct <- IntegrateLayers(object = PTs.combined.sct, method = CCAIntegration, 
                                    orig.reduction = "pca", new.reduction = "integrated.cca",
                                    normalization.method = "SCT", verbose = FALSE)

# Recalcular UMAP/Clusters no espaço agora alinhado ("integrated.cca")
PTs.combined.sct <- RunUMAP(PTs.combined.sct, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.integrated")
PTs.combined.sct <- FindNeighbors(PTs.combined.sct, reduction = "integrated.cca", dims = 1:30)
PTs.combined.sct <- FindClusters(PTs.combined.sct, resolution = 0.9, algorithm = 2)

DimPlot(PTs.combined.sct, reduction = "umap.integrated", label = FALSE, repel = TRUE) + ggtitle("True Integrated Clusters (Clean Data)")
ggsave(paste0(plot_dir, "PT Cluster.png"), width = 12, height = 8, dpi = 300)

# Garantir que a análise diferencial ocorra no assay "SCT"
DefaultAssay(PTs.combined.sct) <- "SCT"
PTs.combined.sct <- PrepSCTFindMarkers(PTs.combined.sct)

PTs_markers <- FindAllMarkers(object = PTs.combined.sct, assay = 'SCT', only.pos = TRUE, logfc.threshold = 0.25)
PTs_top <- PTs_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

sheetlist <- list("allPosMarkers"=PTs_markers,"Top20Markers"=PTs_top)
write_xlsx(sheetlist, paste0(result_dir, 'PTs_markers.xlsx'))

# Validação com DotPlot e Stacked Violin (PTs)
PTs.marker <- c("BRCA1","ESR1","BCL2","PGR","CCNB1","MKI67","AURKA",
                "BAG1", "KRT8","BRIP1", "TUBB", "TUBA1B", "KRT17", "CAV1",
                "IFI27","IGFBP4","TFF1","FTL","FGF13", "SLC3A2")

DotPlot(PTs.combined.sct, features = PTs.marker, group.by = 'seurat_clusters', cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave(paste0(plot_dir, "PT DotPlot.png"), width = 12, height = 8, dpi = 300)

VlnPlot(PTs.combined.sct, features = PTs.marker, stack=T,pt.size=0, flip = T, add.noise = T) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 14), legend.position = 'none')
ggsave(paste0(plot_dir, "PTs_Canonical_Markers_Violin.png"), width = 10, height = 8, dpi = 300)

# Re-anotar os clusters PT com base na avaliação biológica
PTs.combined.sct <- RenameIdents(object = PTs.combined.sct,
                                 "0" = "Células Cancerígenas", "1" = "Células Cancerígenas",
                                 "2" = "Células Cancerígenas", "3" = "Células Cancerígenas",  
                                 "4" = "Células Cancerígenas", "5" = "Células Cancerígenas",  
                                 "6" = "Células Cancerígenas", "7" = "Células Luminais",    
                                 "8" = "Células Cancerígenas", "9" = "Células Cancerígenas",
                                 "10" = "Células Cancerígenas", "11" = "Células Cancerígenas",      
                                 "12" = "Células NK", "13" = "Células Mielóides",
                                 "14" = "Células Epiteliais", "15" = "Células Basais")

PTs.combined.sct$cell.annot <- Idents(PTs.combined.sct)

DimPlot(PTs.combined.sct, reduction = "umap.integrated", label = FALSE, repel = TRUE, label.size = 4)
ggsave(paste0(plot_dir, "PTs_Cluster_anoo.png"), width = 10, height = 8, dpi = 300)

# Heatmap dos 10 principais genes por cluster (Necessita pré-scale no Seurat V5)
PTs_top_10 <- PTs_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
PTs.combined.sct <- ScaleData(PTs.combined.sct, features = unique(PTs_top_10$gene), verbose = FALSE)

DoHeatmap(PTs.combined.sct, features = unique(PTs_top_10$gene), label = FALSE) + theme(text = element_text(face = "bold", size = 8))
ggsave(paste0(plot_dir, "PTs_Heat.png"), width = 10, height = 8, dpi = 300)

saveRDS(PTs.combined.sct, file = paste0(result_dir, "PTs_classified.rds"))


# ------------------------------------------------------------------------------
# 3. TUMORES RECORRENTES (RTs)
# ------------------------------------------------------------------------------
# Separar as 3 amostras de tumor recorrente para avaliação independente
RTs.combined.sct <- subset(TTs.combined.sct, subset = orig.ident %in% c('RT1', 'RT2', 'RT3'))
RTs.combined.sct <- RunPCA(RTs.combined.sct, npcs = 30, verbose = FALSE)

# Integração limpa (Mitigação de batch noise entre os pacientes)
RTs.combined.sct <- IntegrateLayers(object = RTs.combined.sct, method = CCAIntegration, 
                                    orig.reduction = "pca", new.reduction = "integrated.cca",
                                    normalization.method = "SCT", verbose = FALSE)

RTs.combined.sct <- RunUMAP(RTs.combined.sct, reduction = "integrated.cca", dims = 1:14)
RTs.combined.sct <- FindNeighbors(RTs.combined.sct, reduction = "integrated.cca", dims = 1:14)
RTs.combined.sct <- FindClusters(RTs.combined.sct, resolution = 0.9, algorithm = 2)

DimPlot(RTs.combined.sct, label = FALSE, repel = TRUE, label.size = 6)
ggsave(paste0(plot_dir, "RTs_Clusters.png"), width = 10, height = 8, dpi = 300)

# Expressão Diferencial para os Recorrentes
DefaultAssay(RTs.combined.sct) <- "SCT"
RTs.combined.sct <- PrepSCTFindMarkers(RTs.combined.sct)

RTs_markers <- FindAllMarkers(object = RTs.combined.sct, assay = 'SCT', only.pos = TRUE, logfc.threshold = 0.25)
RTs_top <- RTs_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

sheetlist <- list("allPosMarkers" = RTs_markers, "Top20Markers" = RTs_top)
write_xlsx(sheetlist, paste0(result_dir, 'RTs_markers.xlsx'))
saveRDS(RTs.combined.sct, file = paste0(result_dir, 'RTs_clean.rds'))

# Validação Gráfica e Anotação - Marcadores Diagnósticos RT
RTs.marker <- c("MUCL1","KRT19","XBP1","NEAT1","ESR1","CCNB1","MKI67","AURKA",
                "TUBB", "TUBA1B","UBE2T", "PCNA", "MSH6", "AGR2", "PDIA4",
                "KRT80","SLC3A2","CD83","NFKBIA","NFKBIZ")

DotPlot(RTs.combined.sct, features = RTs.marker, group.by = 'seurat_clusters', cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave(paste0(plot_dir, "RTs_DotPlot.png"), width = 10, height = 8, dpi = 300)

VlnPlot(RTs.combined.sct, features = RTs.marker, stack = TRUE, flip = TRUE, pt.size = 0, group.by = "seurat_clusters") + 
  NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(plot_dir, "RTs_violin.png"), width = 10, height = 8, dpi = 300)

# Template de anotação provisório (Mapeie após checar os marcadores)
RTs.combined.sct <- RenameIdents(object = RTs.combined.sct,
                                 "0"  = "Células Cancerígenas", "1"  = "Células Cancerígenas", 
                                 "2"  = "Células Epiteliais", "3"  = "Células Cancerígenas", 
                                 "4"  = "Células Cancerígenas", "5"  = "Células Cancerígenas", 
                                 "6"  = "Células Endoteliais", "7"  = "Células Mielóides", 
                                 "8"  = "Células Cancerígenas", "9"  = "Células Cancerígenas", 
                                 "10" = "Células Cancerígenas", "11" = "Células Cancerígenas",
                                 "12" = "Células Cancerígenas", "13" = "Células Cancerígenas",
                                 "14" = "Células Cancerígenas", "15" = "Células Cancerígenas") 

RTs.combined.sct$cell.annot <- Idents(RTs.combined.sct)

DimPlot(RTs.combined.sct, reduction = "umap", label = FALSE, repel = TRUE, label.size = 4)
ggsave(paste0(plot_dir, "RTs_Cluster_anoo.png"), width = 10, height = 8, dpi = 300)

# Mapa de Calor (Heatmap) Final - Top 8 genes para os RTs
RTs_top_8 <- RTs_top %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
RTs.combined.sct <- ScaleData(RTs.combined.sct, features = unique(RTs_top_8$gene), verbose = FALSE)

DoHeatmap(RTs.combined.sct, features = unique(RTs_top_8$gene), label = FALSE) + theme(text = element_text(face = "bold", size = 8))
ggsave(paste0(plot_dir, "RTs_Heat.png"), width = 10, height = 8, dpi = 300)

saveRDS(RTs.combined.sct, file = paste0(result_dir, "RTs_classified.rds"))