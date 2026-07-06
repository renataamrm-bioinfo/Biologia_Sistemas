# ==============================================================================
# SCRIPT DE VISUALIZAÇÃO: GERAÇÃO DE GRÁFICOS EM PORTUGUÊS
# ==============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(writexl)

# ==============================================================================
# CONFIGURAÇÃO DE DIRETÓRIOS (ATUALIZADO)
# ==============================================================================
base_dir   <- "D:/Onedrive/Renata_Proj/Test_Proj/"

# CORREÇÃO: Apontando diretamente para a sua pasta 'check'
result_dir <- paste0(base_dir, "check/")    
plot_dir   <- paste0(result_dir, "plots_CellTypes/")    # Onde os novos gráficos em PT vão entrar


# Garantir que a pasta de plots em português seja criada se não existir
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# ETAPA 1: Carregar Objetos Processados e Metadados
# ------------------------------------------------------------------------------
print("Carregando objetos RDS...")

# Se você precisar refazer os plots de Pré-QC, precisaremos recalcular os metadados brutos.
# Caso queira focar apenas nos dados limpos e classificados, carregamos abaixo:
NTs.seu <- readRDS(paste0(result_dir, "NTs_classified.rds"))
PTs.combined.sct <- readRDS(paste0(result_dir, "PTs_classified.rds"))
RTs.combined.sct <- readRDS(paste0(result_dir, "RTs_classified.rds"))

# ==============================================================================
# 1. GRÁFICOS DO TECIDO NORMAL (NTs)
# ==============================================================================
print("Gerando gráficos para Tecidos Normais...")

# UMAP Classificado
DimPlot(NTs.seu, label = FALSE, repel = TRUE, label.size = 5) + 
  ggtitle("Tecido Normal Anotado (Pipeline Limpo)") +
  labs(x = "UMAP 1", y = "UMAP 2", fill = "Tipos Celulares") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(paste0(plot_dir, "NT_UMAP_Anotado_PT.png"), width = 10, height = 7, dpi = 300)

# Violino dos Marcadores Canônicos (Traduzindo eixos e títulos)
NTs.marker <- c("KRT15","KRT16","KRT17", "ESRP1","ELF3","RARRES1",
                "TOP2A", "CDK1","MKI67","CENPF", "COL4A6","COL4A5",
                "S100A10", "ID1", "CXCL1", "CXCL8", "CD24")

VlnPlot(NTs.seu, features = NTs.marker, stack = TRUE, pt.size = 0, flip = TRUE, add.noise = TRUE) +
  labs(title = "Expressão de Marcadores Canônicos - Tecido Normal",
       x = "Genes Marcadores", 
       y = "Identidade Celular") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(colour = 'black', size = 12), 
        axis.text.x = element_text(colour = 'black', size = 12, angle = 45, hjust = 1),
        legend.position = 'none')
ggsave(paste0(plot_dir, "NT_Marcadores_Canonicos_Violino_PT.png"), width = 11, height = 8, dpi = 300)


# ==============================================================================
# 2. GRÁFICOS DOS TUMORES PRIMÁRIOS (PTs)
# ==============================================================================
print("Gerando gráficos para Tumores Primários...")

# UMAP Classificado
DimPlot(PTs.combined.sct, reduction = "umap.integrated", label = FALSE, repel = TRUE) + 
  ggtitle("Agrupamentos Celulares Integrados - Tumor Primário (PT)") +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Tipos Celulares") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(paste0(plot_dir, "PT_UMAP_Integrado_PT.png"), width = 10, height = 7, dpi = 300)

# DotPlot de Marcadores Canônicos
PTs.marker <- c("BRCA1","ESR1","BCL2","PGR","CCNB1","MKI67","AURKA",
                "BAG1", "KRT8","BRIP1", "TUBB", "TUBA1B", "KRT17", "CAV1",
                "IFI27","IGFBP4","TFF1","FTL","FGF13", "SLC3A2")

DotPlot(PTs.combined.sct, features = PTs.marker, group.by = 'cell.annot', cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis() +
  labs(title = "Perfil de Expressão Gênica - Tumores Primários",
       x = "Marcadores Moleculares", 
       y = "Grupos Celulares Anotados",
       size = "% Células Expressando",
       color = "Expressão Média") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(paste0(plot_dir, "PT_DotPlot_Marcadores_PT.png"), width = 12, height = 7, dpi = 300)

# Heatmap (Mapa de Calor) dos Top Genes
# Nota: Como re-carregamos o objeto, precisamos garantir que a matriz de escala tenha os genes do heatmap
if(exists("PTs_markers") == FALSE) {
  # Caso não tenha os marcadores no ambiente, lemos direto do Excel gerado no outro script
  PTs_markers <- readxl::read_excel(paste0(result_dir, 'PTs_markers.xlsx'), sheet = "allPosMarkers")
}
PTs_top_10 <- PTs_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
PTs.combined.sct <- ScaleData(PTs.combined.sct, features = unique(PTs_top_10$gene), verbose = FALSE)

DoHeatmap(PTs.combined.sct, features = unique(PTs_top_10$gene), label = FALSE) + 
  ggtitle("Mapa de Calor: Top 10 Genes Marcadores por Cluster (PT)") +
  labs(fill = "Nível de\nExpressão") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        text = element_text(face = "bold", size = 8))
ggsave(paste0(plot_dir, "PT_Heatmap_Top10_PT.png"), width = 11, height = 8, dpi = 300)


# ==============================================================================
# 3. GRÁFICOS DOS TUMORES RECORRENTES (RTs)
# ==============================================================================
print("Gerando gráficos para Tumores Recorrentes...")

# CORREÇÃO: Alterado de "umap.integrated" para "umap" (padrão do seu objeto RT)
DimPlot(RTs.combined.sct, reduction = "umap", label = FALSE, repel = TRUE) + 
  ggtitle("Agrupamentos Celulares Integrados - Tumor Recorrente (RT)") +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Tipos Celulares") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(paste0(plot_dir, "RT_UMAP_Integrado_PT.png"), width = 10, height = 7, dpi = 300)

# DotPlot de Marcadores Diagnósticos RT
RTs.marker <- c("MUCL1","KRT19","XBP1","NEAT1","ESR1","CCNB1","MKI67","AURKA",
                "TUBB", "TUBA1B","UBE2T", "PCNA", "MSH6", "AGR2", "PDIA4",
                "KRT80","SLC3A2","CD83","NFKBIA","NFKBIZ")

# CORREÇÃO: Garantindo que o agrupamento use 'cell.annot' que você salvou
DotPlot(RTs.combined.sct, features = RTs.marker, group.by = 'cell.annot', cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis() +
  labs(title = "Perfil de Expressão Gênica - Tumores Recorrentes",
       x = "Marcadores Moleculares", 
       y = "Grupos Celulares Anotados",
       size = "% Células Expressando",
       color = "Expressão Média") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(paste0(plot_dir, "RT_DotPlot_Marcadores_PT.png"), width = 12, height = 7, dpi = 300)



# ==============================================================================
# NOVA SEÇÃO: PERFIL DE EXPRESSÃO DOS MARCADORES POR CLUSTER CRU (SEURAT_CLUSTERS)
# ==============================================================================
print("Gerando gráficos de marcadores contra os clusters originais (Sem Anotação)...")

# ------------------------------------------------------------------------------
# 1. DotPlot: Marcadores de PT vs Seurat Clusters
# ------------------------------------------------------------------------------
DotPlot(PTs.combined.sct, features = PTs.marker, group.by = 'seurat_clusters', cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis() +
  labs(title = "Expressão de Marcadores de PT por Cluster Original",
       x = "Marcadores Moleculares (PT)", 
       y = "Seurat Clusters (Números)",
       size = "% Células Expressando",
       color = "Expressão Média") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
ggsave(paste0(plot_dir, "PT_DotPlot_Marcadores_vs_Clusters.png"), width = 12, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 2. DotPlot: Marcadores de RT vs Seurat Clusters
# ------------------------------------------------------------------------------
DotPlot(RTs.combined.sct, features = RTs.marker, group.by = 'seurat_clusters', cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis() +
  labs(title = "Expressão de Marcadores de RT por Cluster Original",
       x = "Marcadores Moleculares (RT)", 
       y = "Seurat Clusters (Números)",
       size = "% Células Expressando",
       color = "Expressão Média") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
ggsave(paste0(plot_dir, "RT_DotPlot_Marcadores_vs_Clusters.png"), width = 12, height = 7, dpi = 300)

print("Novos gráficos de marcadores vs clusters salvos com sucesso!")








# ==============================================================================
# GERAÇÃO DE GRÁFICOS DE CLUSTERS (UMAP) - PT E RT ISOLADOS
# ==============================================================================
print("Gerando gráficos de agrupamento (Clusters) para PT e RT...")

# ------------------------------------------------------------------------------
# 1. Gráfico de Clusters para Tumores Primários (PT)
# ------------------------------------------------------------------------------
# Aqui definimos o group.by como 'seurat_clusters' para mostrar os números dos clusters originais
DimPlot(PTs.combined.sct, reduction = "umap.integrated", group.by = "seurat_clusters", label = TRUE, repel = TRUE, label.size = 5) + 
  ggtitle("Agrupamentos Celulares (Clusters) - Tumores Primários (PT)") +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.text = element_text(size = 10))
ggsave(paste0(plot_dir, "Clusters_UMAP_Tumores_Primarios_PT.png"), width = 9, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 2. Gráfico de Clusters para Tumores Recorrentes (RT)
# ------------------------------------------------------------------------------
# Lembrando que no seu RT a redução padrão ficou salva como "umap"
DimPlot(RTs.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, label.size = 5) + 
  ggtitle("Agrupamentos Celulares (Clusters) - Tumores Recorrentes (RT)") +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.text = element_text(size = 10))
ggsave(paste0(plot_dir, "Clusters_UMAP_Tumores_Recorrentes_PT.png"), width = 9, height = 7, dpi = 300)

print("Gráficos de clusters de PT e RT gerados separadamente!")




# ==============================================================================
# SCRIPT DE QC DETALHADO: GRÁFICOS INDIVIDUAIS EM PORTUGUÊS (SEM NTs)
# ==============================================================================
print("Carregando o objeto de tumores totais para métricas de QC...")
TTs.combined.sct <- readRDS(paste0(result_dir, "TTS.rds"))

# 1. Extrair e preparar metadados (Apenas Amostras Tumorais)
metadata_all <- TTs.combined.sct@meta.data %>% 
  filter(!str_detect(orig.ident, "^NT"))

# Função local para cálculo do nMAD (Rigor = 5 conforme seu pipeline)
nMAD_local <- function(x, nmads = 5) {
  xm <- median(x)
  md <- median(abs(x - xm))
  return(ceiling(xm + nmads * md))
}

# Estruturando os cortes para linhas verticais nos gráficos Pré-QC
cortes_ncount <- metadata_all %>%
  group_by(orig.ident) %>%
  summarise(up.cut = nMAD_local(nCount_RNA), down.cut = nMAD_local(nCount_RNA)/20) %>%
  pivot_longer(cols = c(up.cut, down.cut), names_to = "types", values_to = "values")

cortes_nfeature <- metadata_all %>%
  group_by(orig.ident) %>%
  summarise(up.cut = nMAD_local(nFeature_RNA), down.cut = nMAD_local(nFeature_RNA)/20) %>%
  pivot_longer(cols = c(up.cut, down.cut), names_to = "types", values_to = "values")

cortes_mt <- metadata_all %>%
  group_by(orig.ident) %>%
  summarise(values = nMAD_local(percent.mt))

cortes_dbl <- metadata_all %>%
  group_by(orig.ident) %>%
  summarise(values = nMAD_local(scDblFinder.score))

# Criando a tabela de metadados simulando o Pós-QC exato
metadata_final <- metadata_all %>%
  group_by(orig.ident) %>%
  mutate(
    limite_nfeature = nMAD_local(nFeature_RNA),
    limite_ncount = nMAD_local(nCount_RNA),
    limite_mt = min(nMAD_local(percent.mt), 25)
  ) %>%
  ungroup() %>%
  filter(
    nFeature_RNA <= limite_nfeature & nFeature_RNA > limite_nfeature/20,
    nCount_RNA <= limite_ncount & nCount_RNA > limite_ncount/20,
    percent.mt <= limite_mt,
    scDblFinder.class == "singlet"
  )

# ==============================================================================
# BLOCO 1: GRÁFICOS PRÉ-QC (INDIVIDUAIS)
# ==============================================================================
print("Gerando arquivos individuais de Pré-QC...")

# Gráfico 1: Densidade de Contagem de RNA (Pré-QC)
metadata_all %>% 
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
  geom_density(alpha = .4) + scale_x_log10() + theme_classic() +
  geom_vline(data = cortes_ncount, aes(xintercept = values), linetype = "dashed", color = "black") +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "1. Densidade de Contagem de RNA (nCount_RNA) - Pré-QC",
       x = "Contagem de RNA (log10)", y = "Densidade", fill = "Amostra", color = "Amostra") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_01_Pre_Contagem_Densidade.png"), width = 10, height = 5, dpi = 300)

# Gráfico 2: Densidade de Genes Detectados (Pré-QC)
metadata_all %>% 
  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
  geom_density(alpha = .4) + scale_x_log10() + theme_classic() +
  geom_vline(data = cortes_nfeature, aes(xintercept = values), linetype = "dashed", color = "black") +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "2. Densidade de Genes Detectados (nFeature_RNA) - Pré-QC",
       x = "Número de Genes (log10)", y = "Densidade", fill = "Amostra", color = "Amostra") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_02_Pre_Genes_Densidade.png"), width = 10, height = 5, dpi = 300)

# Gráfico 3: Densidade de Porcentagem Mitocondrial (Pré-QC)
metadata_all %>% 
  ggplot(aes(color = orig.ident, x = percent.mt, fill = orig.ident)) + 
  geom_density(alpha = .4) + theme_classic() +
  geom_vline(data = cortes_mt, aes(xintercept = values), linetype = "dashed", color = "black") +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "3. Densidade de Porcentagem Mitocondrial (percent.mt) - Pré-QC",
       x = "% Mitocondrial (Linear)", y = "Densidade", fill = "Amostra", color = "Amostra") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_03_Pre_Mitocondrial_Densidade.png"), width = 10, height = 5, dpi = 300)

# Gráfico 4: Densidade do Score de Doublet (Pré-QC)
metadata_all %>% 
  ggplot(aes(color = orig.ident, x = scDblFinder.score, fill = orig.ident)) + 
  geom_density(alpha = .4) + theme_classic() +
  geom_vline(data = cortes_dbl, aes(xintercept = values), linetype = "dashed", color = "black") +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "4. Densidade do Score de Doublet (scDblFinder) - Pré-QC",
       x = "Score de Doublet", y = "Densidade", fill = "Amostra", color = "Amostra") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_04_Pre_Doublet_Densidade.png"), width = 10, height = 5, dpi = 300)

# Gráfico 5: Dispersão nCount vs nFeature (Pré-QC)
metadata_all %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point(size = 0.6, alpha = 0.6) + scale_colour_gradient(low = "gray85", high = "firebrick") +
  geom_smooth(se = TRUE, level = 0.9, color = "blue", method = "gam") + scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(data = cortes_ncount, aes(xintercept = values), linetype = "dashed", color = "black") +
  geom_hline(data = cortes_nfeature, aes(yintercept = values), linetype = "dashed", color = "black") +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "5. Dispersão de nCount_RNA vs nFeature_RNA - Pré-QC",
       x = "Contagem de RNA (log10)", y = "Número de Genes (log10)", color = "% Mitocondrial") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_05_Pre_Dispersao_Correlação.png"), width = 11, height = 6, dpi = 300)


# ==============================================================================
# BLOCO 2: GRÁFICOS PÓS-QC (INDIVIDUAIS)
# ==============================================================================
print("Gerando arquivos individuais de Pós-QC...")

# Gráfico 6: Densidade de Contagem de RNA (Pós-QC)
metadata_final %>% 
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
  geom_density(alpha = .4) + scale_x_log10() + theme_classic() +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "6. Densidade de Contagem de RNA (nCount_RNA) - Pós-QC",
       x = "Contagem de RNA (log10)", y = "Densidade", fill = "Amostra", color = "Amostra") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_06_Pos_Contagem_Densidade.png"), width = 10, height = 5, dpi = 300)

# Gráfico 7: Densidade de Genes Detectados (Pós-QC)
metadata_final %>% 
  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
  geom_density(alpha = .4) + scale_x_log10() + theme_classic() +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "7. Densidade de Genes Detectados (nFeature_RNA) - Pós-QC",
       x = "Número de Genes (log10)", y = "Densidade", fill = "Amostra", color = "Amostra") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_07_Pos_Genes_Densidade.png"), width = 10, height = 5, dpi = 300)

# Gráfico 8: NOVO! Dispersão nCount vs nFeature (Pós-QC)
metadata_final %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point(size = 0.6, alpha = 0.6) + scale_colour_gradient(low = "gray85", high = "firebrick") +
  geom_smooth(se = TRUE, level = 0.9, color = "blue", method = "gam") + scale_x_log10() + scale_y_log10() + theme_classic() +
  facet_wrap(~orig.ident, ncol = 3) +
  labs(title = "8. Dispersão de nCount_RNA vs nFeature_RNA - Pós-QC",
       x = "Contagem de RNA (log10)", y = "Número de Genes (log10)", color = "% Mitocondrial") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
ggsave(paste0(plot_dir, "QC_08_Pos_Dispersao_Correlacao.png"), width = 11, height = 6, dpi = 300)

# Gráfico 9: Comparação da Quantidade de Células (Pré vs Pós)
df_contagem <- rbind(
  data.frame(Amostra = metadata_all$orig.ident, Estado = "Pré-QC"),
  data.frame(Amostra = metadata_final$orig.ident, Estado = "Pós-QC")
)

df_contagem %>%
  ggplot(aes(x = Amostra, fill = Estado)) +
  geom_bar(position = "dodge") +
  theme_classic() +
  labs(title = "9. Impacto do Controle de Qualidade no Número Total de Células",
       x = "Amostra Tumoral", y = "Número de Células", fill = "Etapa do Pipeline") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(plot_dir, "QC_09_Comparativo_Quantidade_Celulas.png"), width = 8, height = 5, dpi = 300)

print("Todos os gráficos individuais foram salvos com sucesso!")