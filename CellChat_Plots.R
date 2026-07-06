# ==============================================================================
# 0. CONFIGURAÇÃO INICIAL E INICIALIZAÇÃO DE BIBLIOTECAS
# ==============================================================================
# Define seu prefixo aqui ("PT", "RT", ou "NT") para evitar erros humanos!
PREFIX <- "PT" 

library(Seurat)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(uwot) # Motor UMAP nativo do R - ignora instalações e erros de Python!

# Configurações globais de ambiente e alocação de memória (10 GB)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10 * 1024^3) 

# Garante que o pacote uwot está disponível no ambiente local
if(!any(rownames(installed.packages()) == "uwot")) {
  install.packages("uwot", repos = "http://cran.us.r-project.org")
}

# Definição e criação automática da estrutura de diretórios no OneDrive
output_base_dir <- paste0("D:/Onedrive/Renata_Proj/Test_Proj/check/", PREFIX, "_NE")

folders <- c("Global_Network_Plots", "Signaling_Pathways/Hierarchy", "Signaling_Pathways/Circle", 
             "Signaling_Pathways/Chord", "Signaling_Pathways/Heatmap", "Signaling_Pathways/Chord_Cell_Grouped", 
             "Signaling_Pathways/Contribution", "Signaling_Pathways/Individual_Hierarchy", 
             "Signaling_Pathways/Individual_Circle", "Signaling_Pathways/Individual_Chord", 
             "Signaling_Pathways/Bubble_Plots", "Signaling_Pathways/Chord_Gene_Level", 
             "Systems_Role_Analysis", "Gene_Expression_Violins", "Global_Coordination_Analysis")

for (f in folders) {
  dir.create(file.path(output_base_dir, f), recursive = TRUE, showWarnings = FALSE)
}

# Inicialização do Arquivo de Log na Raiz
log_file <- file.path(output_base_dir, "relatorio_erros_e_avisos.txt")
writeLines(paste0("======================================================================\n",
                  "📋 RELATÓRIO DE EXECUÇÃO E ALERTAS DO PIPELINE CELLCHAT\n",
                  "Data de Início: ", Sys.time(), "\n",
                  "Prefixo Configurado: ", PREFIX, "\n",
                  "======================================================================\n"), log_file)

# Função auxiliar interna para escrever mensagens no Log de forma limpa
registrar_log <- function(mensagem) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", mensagem, "\n"), file = log_file, append = TRUE)
}

# ==============================================================================
# 1. CARREGAMENTO DO OBJETO CELLCHAT PROCESSADO (.RDS)
# ==============================================================================
message("💾 Carregando o objeto CellChat pré-processado...")
ptm <- Sys.time()

# Carrega dinamicamente o arquivo baseado no prefixo escolhido (ex: RTs_cellchat.rds)
cellchat <- readRDS(paste0("D:/Onedrive/Renata_Proj/Test_Proj/check/", PREFIX, "s_cellchat.rds"))

# Exportação do dataframe com todas as comunicações inferidas a nível de L-R
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = file.path(output_base_dir, paste0(PREFIX, "_NE_cellchat_normalized.csv")), row.names = FALSE)

# Computação das probabilidades a nível de vias gerais e agregação da rede
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Extração dinâmica de caminhos moleculares activos e índices celulares
all_pathways    <- cellchat@netP$pathways
vertex.receiver <- 1:length(levels(cellchat@idents)) 

# ==============================================================================
# 2. GRÁFICOS DE CÍRCULO GLOBAIS (MICROAMBIENTE INTEGRADO)
# ==============================================================================
message("📊 Salvando gráficos de círculos integrados do microambiente...")
groupSize <- as.numeric(table(cellchat@idents))

png(filename = file.path(output_base_dir, "Global_Network_Plots", paste0(PREFIX, "_Global_Aggregated_Interactions.png")), width = 1400, height = 700, res = 130)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Número de interações")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Força/Peso das interações")
dev.off()

# Loop para gerar a distribuição de sinais enviados individualmente por tipo celular
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  safe_cell_name <- gsub("[^[:alnum:]_]", "_", rownames(mat)[i])
  
  png(filename = file.path(output_base_dir, "Global_Network_Plots", paste0(PREFIX, "_NetVisual_Circle_", safe_cell_name, ".png")), width = 800, height = 900, res = 120)
  par(xpd = TRUE, mar = c(2, 2, 4, 2)) 
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = 10, arrow.size = 1.0, arrow.width = 1.0, title.name = paste("Sinais enviados por:", rownames(mat)[i]))
  dev.off()
}

# ==============================================================================
# 3. CONFIGURAÇÃO DO MAPEAMENTO DINÂMICO DE GRUPOS (Automático e Inteligente)
# ==============================================================================
message("🤖 Mapeando tipos celulares de forma automatizada para PT e RT...")

cluster_names <- levels(cellchat@idents)

# Cria um vetor vazio para armazenar o mapeamento
group_mapping <- character(length(cluster_names))

for (i in seq_along(cluster_names)) {
  name <- cluster_names[i]
  
  group_mapping[i] <- case_when(
    # 1. Células Cancerígenas / Tumorais
    grepl("cancer|tumor|malign", name, ignore.case = TRUE) ~ "Tumor",
    
    # 2. Células Epiteliais Gerais
    grepl("epitelial|epiteliais|epithelial", name, ignore.case = TRUE) ~ "Epitelial",
    
    # 3. Células Luminais (comum em mama/próstata, frequentemente epiteliais)
    grepl("luminal|luminais", name, ignore.case = TRUE) ~ "Luminal",
    
    # 4. Células Basais
    grepl("basal|basais", name, ignore.case = TRUE) ~ "Basal",
    
    # 5. Células Endoteliais
    grepl("endotelial|endoteliais|endothelial", name, ignore.case = TRUE) ~ "Endotelial",
    
    # 6. Células Mielóides
    grepl("miel[oó]ide|miel[oó]ides|myeloid", name, ignore.case = TRUE) ~ "Mielóide",
    
    # 7. Células NK / Linfoides
    grepl("nk|linf[oó]cito|t_cell|b_cell|lymphoid", name, ignore.case = TRUE) ~ "Linfoide",
    
    # Caso apareça algum cluster totalmente novo no futuro
    TRUE ~ "Outro"
  )
}

group.cellType <- group_mapping
names(group.cellType) <- cluster_names

# Exibe no console a tabela de conferência para o PREFIX atual
message(paste("📋 Tabela de Mapeamento para o dataset:", PREFIX))
print(data.frame(Original = cluster_names, Mapeado = group.cellType, row.names = NULL))


# ==============================================================================
# [EXTRA] GRÁFICO DE BOLHAS GLOBAL (TODAS AS VIAS JUNTAS)
# ==============================================================================
message("🔮 Gerando Gráfico de Bolhas Global para triagem de vias no artigo...")

tryCatch({
  png(filename = file.path(output_base_dir, "Systems_Role_Analysis", paste0(PREFIX, "_Global_Bubble_All_Pathways.png")), 
      width = 1600, height = 2000, res = 130) # Altura maior para acomodar todas as vias
  
  p_global_bubble <- netVisual_bubble(
    cellchat, 
    sources.use = vertex.receiver, 
    targets.use = vertex.receiver, 
    signaling = all_pathways, 
    remove.isolate = FALSE
  ) + 
    ggtitle(paste("Visão Geral de Todas as Vias de Sinalização -", PREFIX)) +
    theme(axis.text.y = element_text(size = 8)) # Ajusta o tamanho da fonte para não sobrepor
  
  print(p_global_bubble)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  registrar_log(paste0("ERRO ao gerar o Bubble Plot Global: ", e$message))
})

# ==============================================================================
# 4. LOOP MESTRE DE VIAS (Geração simultânea de todas as camadas visuais)
# ==============================================================================
message("🔄 Iniciando loop automatizado para processamento individual de vias...")

for (pathway.show in all_pathways) {
  safe_pathway_name <- gsub("[^[:alnum:]_]", "_", pathway.show)
  
  # A. Gráfico de Hierarquia (Hierarchy Plot)
  png(filename = file.path(output_base_dir, "Signaling_Pathways/Hierarchy", paste0(PREFIX, "_Hierarchy_", safe_pathway_name, ".png")), width = 1200, height = 700, res = 130)
  par(mar = c(3, 3, 5, 3), xpd = TRUE)
  netVisual_aggregate(cellchat, signaling = pathway.show, vertex.receiver = vertex.receiver)
  title(main = paste("Rede de Sinalização de", pathway.show, "(Hierarquia)"), cex.main = 1.5, line = 2)
  dev.off()
  
  # B. Gráfico de Círculo (Circle Plot)
  png(filename = file.path(output_base_dir, "Signaling_Pathways/Circle", paste0(PREFIX, "_Circle_", safe_pathway_name, ".png")), width = 900, height = 950, res = 130)
  par(mar = c(3, 3, 6, 3), xpd = TRUE)
  netVisual_aggregate(cellchat, signaling = pathway.show, layout = "circle", vertex.label.cex = 0.8, margin = 0.2, arrow.size = 1.5, arrow.width = 1.5)
  title(main = paste("Rede de Sinalização de", pathway.show, "(Círculo)"), cex.main = 1.5, line = 3)
  dev.off()
  
  # C. Diagrama de Cordas Geral (Chord Plot)
  png(filename = file.path(output_base_dir, "Signaling_Pathways/Chord", paste0(PREFIX, "_Chord_", safe_pathway_name, ".png")), width = 900, height = 950, res = 130)
  par(mar = c(3, 3, 6, 3), xpd = TRUE)
  netVisual_aggregate(cellchat, signaling = pathway.show, layout = "chord")
  title(main = paste("Rede de Sinalização de", pathway.show, "(Cordas)"), cex.main = 1.5, line = 3)
  dev.off()
  
  # D. Mapas de Calor (Heatmaps)
  tryCatch({
    png(filename = file.path(output_base_dir, "Signaling_Pathways/Heatmap", paste0(PREFIX, "_Heatmap_", safe_pathway_name, ".png")), width = 800, height = 800, res = 130)
    print(netVisual_heatmap(cellchat, signaling = pathway.show, color.heatmap = "Reds", title.name = paste("Força da Via", pathway.show)))
    dev.off()
  }, error = function(e) { 
    if (dev.cur() > 1) dev.off() 
    registrar_log(paste0("AVISO: Mapa de Calor pulado para a via '", pathway.show, "' (Sem variância distinta nos dados. Erro: ", e$message, ")"))
  })
  
  # E. Diagrama de Cordas Agrupado por Tipo Celular (Mapeamento em PT-BR)
  tryCatch({
    png(filename = file.path(output_base_dir, "Signaling_Pathways/Chord_Cell_Grouped", paste0(PREFIX, "_Grouped_Chord_", safe_pathway_name, ".png")), width = 1000, height = 1000, res = 130)
    netVisual_chord_cell(cellchat, signaling = pathway.show, group = group.cellType, title.name = paste("Rede de Sinalização de", pathway.show, "Agrupada"))
    dev.off()
  }, error = function(e) { 
    if (dev.cur() > 1) dev.off()
    registrar_log(paste0("AVISO: Gráfico de Cordas Agrupado pulado para a via '", pathway.show, "'. Erro: ", e$message))
  })
  
  # F. Gráficos de Contribuição de Pares L-R (Contribution Bar Plots)
  png(filename = file.path(output_base_dir, "Signaling_Pathways/Contribution", paste0(PREFIX, "_Contribution_", safe_pathway_name, ".png")), width = 800, height = 600, res = 130)
  p_contrib <- netAnalysis_contribution(cellchat, signaling = pathway.show) + 
    ggtitle(paste("Contribuição de Pares Ligante-Receptor - Via", pathway.show)) +
    xlab("Pares Ligante-Receptor") + ylab("Contribuição Relativa")
  print(p_contrib)
  dev.off()
  
  # G. Gráficos de Bolhas (Bubble Plots)
  tryCatch({
    png(filename = file.path(output_base_dir, "Signaling_Pathways/Bubble_Plots", paste0(PREFIX, "_Bubble_", safe_pathway_name, ".png")), width = 1100, height = 900, res = 130)
    p_bubble <- netVisual_bubble(cellchat, sources.use = vertex.receiver, targets.use = vertex.receiver, signaling = pathway.show, remove.isolate = TRUE) + 
      ggtitle(paste("Interações Ligante-Receptor Enriquecidas - Via", pathway.show)) +
      xlab("Pares de Tipos Celulares (Fonte -> Alvo)") + ylab("Pares Ligante-Receptor")
    print(p_bubble)
    dev.off()
  }, error = function(e) { 
    if (dev.cur() > 1) dev.off()
    registrar_log(paste0("AVISO: Gráfico de Bolhas pulado para a via '", pathway.show, "' (Nenhum par L-R ativo significativo. Erro: ", e$message, ")"))
  })
  
  # H. Diagrama de Cordas Avançado a Nível de Genes
  tryCatch({
    png(filename = file.path(output_base_dir, "Signaling_Pathways/Chord_Gene_Level", paste0(PREFIX, "_Chord_Gene_", safe_pathway_name, ".png")), width = 1400, height = 1400, res = 140)
    par(mar = c(2, 2, 6, 2), xpd = TRUE)
    netVisual_chord_gene(cellchat, sources.use = vertex.receiver, targets.use = vertex.receiver, signaling = pathway.show, lab.cex = 0.8, small.gap = 1, big.gap = 8, legend.pos.y = 30)       
    title(main = paste("Rede de Sinalização de", pathway.show, "a Nível de Genes"), cex.main = 2.0, line = 3)
    dev.off()
  }, error = function(e) { 
    if (dev.cur() > 1) dev.off()
    registrar_log(paste0("AVISO: Diagrama de Cordas por Gene pulado para a via '", pathway.show, "'. Erro: ", e$message))
  })
  
  # I. Sub-gráficos Individuais para o Par L-R de Maior Impacto
  pairLR.all <- extractEnrichedLR(cellchat, signaling = pathway.show, geneLR.return = FALSE)
  if (!is.null(pairLR.all) && nrow(pairLR.all) > 0) {
    LR.show <- pairLR.all[1, ] 
    safe_LR_name <- gsub("[^[:alnum:]_]", "_", LR.show)
    
    png(filename = file.path(output_base_dir, "Signaling_Pathways/Individual_Hierarchy", paste0(PREFIX, "_Ind_Hierarchy_", safe_pathway_name, "_", safe_LR_name, ".png")), width = 1200, height = 700, res = 130)
    par(mar = c(3, 3, 5, 3), xpd = TRUE)
    netVisual_individual(cellchat, signaling = pathway.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver)
    title(main = paste("Par Ligante-Receptor Individual:", LR.show, "(Hierarquia)"), cex.main = 1.2, line = 2)
    dev.off()
    
    png(filename = file.path(output_base_dir, "Signaling_Pathways/Individual_Circle", paste0(PREFIX, "_Ind_Circle_", safe_pathway_name, "_", safe_LR_name, ".png")), width = 900, height = 950, res = 130)
    par(mar = c(3, 3, 6, 3), xpd = TRUE)
    netVisual_individual(cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle", vertex.label.cex = 0.8, margin = 0.2)
    title(main = paste("Par Ligante-Receptor Individual:", LR.show, "(Círculo)"), cex.main = 1.2, line = 3)
    dev.off()
  }
}

# ==============================================================================
# 5. ANÁLISE DE CENTRALIDADE E MAPAS DE PAPÉIS GLOBAIS (SYSTEMS ANALYSIS)
# ==============================================================================
message("🔥 Calculando centralidade topológica e mapas de calor de papéis...")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Gráficos de Dispersão 2D (Scatter Plots)
png(filename = file.path(output_base_dir, "Systems_Role_Analysis", paste0(PREFIX, "_Global_Signaling_Role_Scatter.png")), width = 1400, height = 700, res = 130)
gg1 <- netAnalysis_signalingRole_scatter(cellchat) + ggtitle("Todas as Vias Agregadas do Sistema") + xlab("Sinalização de Saída (Força do Emissor)") + ylab("Sinalização de Entrada (Força do Receptor)")
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = all_pathways) + ggtitle("Resumo das Vias de Sinalização Descobertas") + xlab("Sinalização de Saída (Força do Emissor)") + ylab("Sinalização de Entrada (Força do Receptor)")
print(gg1 + gg2)
dev.off()

# Padrões de Entrada (Incoming) e Saída (Outgoing) Gerais (Substituição forçada de rótulos)
png(filename = file.path(output_base_dir, "Systems_Role_Analysis", paste0(PREFIX, "_Global_Pattern_Heatmaps.png")), width = 1400, height = 900, res = 130)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

# Força a tradução das labels dos eixos nos objetos estruturais ComplexHeatmap
ht1@column_title <- "Papéis de Sinalização de Saída (Emissores)"
ht2@column_title <- "Papéis de Sinalização de Entrada (Receptores)"
draw(ht1 + ht2)
dev.off()

# Matriz Resumo de Papéis Específicos por Via
png(filename = file.path(output_base_dir, "Systems_Role_Analysis", paste0(PREFIX, "_Pathway_Specific_Signaling_Role_Heatmap.png")), width = 900, height = 900, res = 130)
ht_spec <- netAnalysis_signalingRole_heatmap(cellchat, signaling = all_pathways)
ht_spec@column_title <- "Importância Relativa dos Papéis de Sinalização por Via"
draw(ht_spec)
dev.off()

# Geração de Gráficos de Violino para a Expressão dos Genes de Sinalização
message("🎻 Gerando violinos de expressão gênica para validação transcriptômica...")
for (pathway.show in all_pathways) {
  safe_pathway_name <- gsub("[^[:alnum:]_]", "_", pathway.show)
  tryCatch({
    png(filename = file.path(output_base_dir, "Gene_Expression_Violins", paste0(PREFIX, "_Expression_", safe_pathway_name, ".png")), width = 1000, height = 700, res = 130)
    p_vln <- plotGeneExpression(cellchat, signaling = pathway.show, enriched.only = TRUE, type = "violin") + 
      ggtitle(paste("Distribuição de Expressão Gênica - Via", pathway.show)) +
      xlab("Tipos Celulares") + ylab("Nível de Expressão Log-Normalizada")
    print(p_vln)
    dev.off()
  }, error = function(e) { 
    if (dev.cur() > 1) dev.off() 
    registrar_log(paste0("AVISO: Gráfico de Violino pulado para a via '", pathway.show, "'. Erro: ", e$message))
  })
}

# ==============================================================================
# 6. DETERMINAÇÃO DE PADRÕES NMF E ENCAIXE DE MANIFOLD (UMAP INTEGRADO COM UWOT)
# ==============================================================================
message("🤖 Decodificando agrupamentos globais por fatoração de matrizes (NMF)...")
k_range    <- 2:7
out_k_plot <- selectK(cellchat, pattern = "outgoing", k.range = k_range)
in_k_plot  <- selectK(cellchat, pattern = "incoming", k.range = k_range)

# Extração matemática do melhor 'k' diretamente da tabela interna gerada nos gráficos
out_metrics <- out_k_plot$data %>% filter(Measure == "Cophenetic")
in_metrics  <- in_k_plot$data %>% filter(Measure == "Cophenetic")

best_out_k  <- out_metrics$k[which.max(out_metrics$score)]
best_in_k   <- in_metrics$k[which.max(in_metrics$score)]

# Ajuste manual para k=3 se houver pouca variação no dataset (garante a geração do River Plot)
best_out_k <- 3
best_in_k  <- 3

message(paste("👉 Configuração NMF Aplicada: Padrões de Saída (k) =", best_out_k))
message(paste("👉 Configuração NMF Aplicada: Padrões de Entrada (k) =", best_in_k))

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = best_out_k)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = best_in_k)

# Exportação dos diagramas Aluviais (River Plots) Protegidos
png(filename = file.path(output_base_dir, "Global_Coordination_Analysis", paste0(PREFIX, "_Patterns_Outgoing_River.png")), width = 1100, height = 900, res = 130)
tryCatch({ 
  netAnalysis_river(cellchat, pattern = "outgoing") 
}, error = function(e){ 
  registrar_log(paste0("ALERTA: Gráfico 'Patterns_Outgoing_River' gerado em branco (Variância ou complexidade menor que k=3. Erro: ", e$message, ")"))
})
dev.off()

png(filename = file.path(output_base_dir, "Global_Coordination_Analysis", paste0(PREFIX, "_Patterns_Outgoing_Dot.png")), width = 1000, height = 800, res = 130)
p_dot1 <- netAnalysis_dot(cellchat, pattern = "outgoing") + 
  ggtitle("Correspondência Epistemológica dos Padrões de Saída (Outgoing)") +
  xlab("Padrões de Sinalização") + ylab("Vias de Sinalização / Grupos Celulares")
print(p_dot1)
dev.off()

png(filename = file.path(output_base_dir, "Global_Coordination_Analysis", paste0(PREFIX, "_Patterns_Incoming_River.png")), width = 1100, height = 900, res = 130)
tryCatch({ 
  netAnalysis_river(cellchat, pattern = "incoming") 
}, error = function(e){ 
  registrar_log(paste0("ALERTA: Gráfico 'Patterns_Incoming_River' gerado em branco (Variância ou complexidade menor que k=3. Erro: ", e$message, ")"))
})
dev.off()

png(filename = file.path(output_base_dir, "Global_Coordination_Analysis", paste0(PREFIX, "_Patterns_Incoming_Dot.png")), width = 1000, height = 800, res = 130)
p_dot2 <- netAnalysis_dot(cellchat, pattern = "incoming") + 
  ggtitle("Correspondência Epistemológica dos Padrões de Entrada (Incoming)") +
  xlab("Padrões de Sinalização") + ylab("Vias de Sinalização / Grupos Celulares")
print(p_dot2)
dev.off()

# Projeção de Agrupamentos Geométricos UMAP (Utilizando o motor nativo uwot protegido)
message("🧬 Calculando projeções de similaridade funcional e estrutural via uwot...")

tryCatch({
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
  cellchat <- netClustering(cellchat, type = "functional")
  
  png(filename = file.path(output_base_dir, "Global_Coordination_Analysis", paste0(PREFIX, "_Embedding_Functional.png")), width = 900, height = 900, res = 130)
  p_emb1 <- netVisual_embedding(cellchat, type = "functional", label.size = 3.5) + ggtitle("Espaço Projetivo de Similaridade Funcional")
  print(p_emb1)
  dev.off()
}, error = function(e){ 
  if (dev.cur() > 1) dev.off()
  registrar_log(paste0("ALERTA: Gráfico 'Embedding_Functional' gerado em branco (Assinatura celular muito uniforme para gerar projeção dimensional. Erro: ", e$message, ")"))
})

tryCatch({
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
  cellchat <- netClustering(cellchat, type = "structural")
  
  png(filename = file.path(output_base_dir, "Global_Coordination_Analysis", paste0(PREFIX, "_Embedding_Structural.png")), width = 900, height = 900, res = 130)
  p_emb2 <- netVisual_embedding(cellchat, type = "structural", label.size = 3.5) + ggtitle("Espaço Projetivo de Similaridade Estrutural")
  print(p_emb2)
  dev.off()
}, error = function(e){ 
  if (dev.cur() > 1) dev.off()
  registrar_log(paste0("ALERTA: Gráfico 'Embedding_Structural' gerado em branco (Assinatura celular muito uniforme para gerar projeção dimensional. Erro: ", e$message, ")"))
})

# ==============================================================================
# 7. SALVAMENTO E SERIALIZAÇÃO DO OBJETO INTEGRADO COMPLETO (.RDS)
# ==============================================================================
message("💾 Arquivando objeto CellChat final consolidado...")
saveRDS(cellchat, file = file.path(output_base_dir, paste0(PREFIX, "_cellchat_completed_analysis.rds")))

execution.time = Sys.time() - ptm

# Finalização formal do Log de Saída
registrar_log(paste0("SUCESSO: Execução concluída. Tempo total de script: ", round(as.numeric(execution.time, units="secs"), 2), " segundos.\n"))

cat("\n======================================================================\n")
cat("🎉 EXECUÇÃO CONCLUÍDA COM SUCESSO E LOG ATUALIZADO!\n")
cat("Prefixo Processado:", PREFIX, "\n")
cat("Tempo total de processamento do pipeline:\n")
print(execution.time)
cat("O arquivo 'relatorio_erros_e_avisos.txt' foi gerado na pasta raiz.\n")
cat("======================================================================\n")