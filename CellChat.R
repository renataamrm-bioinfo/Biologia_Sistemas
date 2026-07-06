# Carregar bibliotecas necessárias
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Usar ambiente virtual do Python, se necessário
reticulate::use_python("C:/Users/clayt/OneDrive/Documentos/.virtualenvs/r-reticulate/Scripts", required = TRUE)

############# FAZER LOAD

data.dir <- 'C:/Users/clayt/OneDrive/Documentos/CellCommunication/ER+/NE'
setwd(data.dir)

# Carregar os objetos salvos no seu ambiente do R
# Isso irá carregar os objetos com os nomes que eles tinham quando foram salvos: 'cellchat' e 'object.list'
load("C:/Users/clayt/OneDrive/Documentos/CellCommunication/ER+/NE/cellchat_BRCA_Pembro_On_ER_NE.rds")
cellchat<- readRDS("C:/Users/clayt/OneDrive/Documentos/CellCommunication/ER+/NE/cellchat_BRCA_Pembro_On_ER_NE.rds")

###############################################################################
# Medir tempo de execução
ptm = Sys.time()

# Carregar matriz de contagem scRNA-seq
data.input.raw <- readRDS("C:/Users/clayt/OneDrive/Documentos/R_workstation/scRNAseq-analysis/1863-counts_cells_cohort1.rds")

# Carregar metadados
meta <- read.csv("C:/Users/clayt/OneDrive/Documentos/R_workstation/scRNAseq-analysis/1872-BIOKEY_metaData_cohort1_web.csv", row.names = 1)

# Adicionar coluna de labels com tipos celulares
meta$labels <- meta$cellType

# Verificar quais células pertencem ao grupo "treatment_naive"
cell.use <- rownames(meta)[meta$timepoint == "On" & meta$expansion == "NE" & meta$BC_type == "TNBC"]

# Subset dos dados
data.input.raw <- data.input.raw[, cell.use]
meta <- meta[cell.use, ]

# Criar objeto Seurat para normalização
seurat.obj <- CreateSeuratObject(counts = data.input.raw, meta.data = meta)

# Normalizar os dados (log-normalização)
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# (Opcional) Escalar os dados e encontrar variáveis  - Apenas para análises no Seurat
seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- FindVariableFeatures(seurat.obj)

# Extrair matriz normalizada
data.input <- GetAssayData(seurat.obj, slot = "data")  # log-normalized matrix

# (opcional) meta <- seurat.obj@meta.data  # já foi usado acima, mas pode atualizar aqui também

# Verificar tipos celulares presentes
unique(meta$labels)

# meta$samples <- meta$patient_id  # ou outro identificador relevante
meta$samples <- NULL

# Criar objeto CellChat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  APOE+ FIB, FBN1+ FIB, COL11A1+ FIB, Inflam. FIB, cDC1, cDC2, LC, Inflam. DC, TC, Inflam. TC, CD40LG+ TC, NKT

#Set the ligand-receptor interaction database
#When analyzing human samples, use the database CellChatDB.human; when analyzing mouse samples, use the database CellChatDB.mouse. CellChatDB categorizes ligand-receptor pairs into different types, including “Secreted Signaling”, “ECM-Receptor”, “Cell-Cell Contact” and “Non-protein Signaling”. By default, the “Non-protein Signaling” are not used.

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

### FORMAS DE APLICAÇÃO:

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use

### Preprocessing the expression data for cell-cell communication analysis:

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
library(future)
plan("sequential")  # Executa tudo em um único processo

options(future.globals.maxSize = 10 * 1024^3)  # 10 GB

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 13.20763
###não consegui fazer abaixo:
#Esse passo não é obrigatório, mas pode aumentar a robustez biológica da análise ao restringir as inferências a interações com suporte de PPI.
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#cellchat <- projectData(cellchat, PPI.human)

### PARTE 2 - inferência da rede de comunicação célula-célula

#Computar a probabilidade de comunicação e inferir rede de comunicação célula-célula

ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-02-14 00:32:35.767285]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-02-14 00:33:13.121225]"

###FILTRAGEM PARA MÍNIMO DE 10 CÉLULAS::
#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)

####EXTRAIR COMO DATAFRAME:
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

##salvar o df.net
output_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE"

# 1. Corrected file.path syntax (no leading slash)
csv_path <- file.path(output_dir, "RT_NE_cellchat_normalized.csv")

# 2. Force R to create the directory if it doesn't exist/is hidden by OneDrive
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 3. Try writing the CSV again
write.csv(df.net, file = csv_path, row.names = FALSE)
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

### Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.

cellchat <- computeCommunProbPathway(cellchat)

###Calculate the aggregated cell-cell communication network:
#CellChat calculates the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. Users can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.

cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#CellChat can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also 
#control the parameter edge.weight.max so that we can compare edge weights between differet networks.

# 1. Define and prep your output directory
output_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Extract your data
mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))

# 3. The Loop
for (i in 1:nrow(mat)) {
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  cell_name <- rownames(mat)[i]
  safe_cell_name <- gsub("[^[:alnum:]_]", "_", cell_name)
  file_name <- paste0("NetVisual_Circle_", safe_cell_name, ".png")
  full_path <- file.path(output_dir, file_name)
  
  png(filename = full_path, width = 800, height = 900, res = 120)
  par(xpd = TRUE, mar = c(2, 2, 4, 2)) 
  
  # Draw the plot with boosted arrow and edge scaling
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = TRUE, 
                   edge.weight.max = 10,  # FIX 1: Max width of the edges (default is usually smaller)
                   arrow.size = 1.0,      # FIX 2: Increases the actual arrow head size (default is usually 0.5 - 1)
                   arrow.width = 1.0,     # FIX 3: Makes the arrow head wider
                   title.name = cell_name)
  
  dev.off()
}





### PARTE 3: Visualização da rede de comunicação célula-célula

#Here we take input of one signaling pathway as an example. All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways


# 1. Define paths and get all available pathways
output_base_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE/Signaling_Pathways"
all_pathways <- cellchat@netP$pathways

# 2. Define your target receiver cells for the hierarchy plot
vertex.receiver <- seq(1, 4) 

# 3. Create subdirectories for each layout type to keep files organized
modes <- c("Hierarchy", "Circle", "Chord", "Heatmap")
for (mode in modes) {
  dir.create(file.path(output_base_dir, mode), recursive = TRUE, showWarnings = FALSE)
}

# 4. Master Loop through every pathway
for (pathway.show in all_pathways) {
  
  # Create a Windows-safe file name for the current pathway
  safe_pathway_name <- gsub("[^[:alnum:]_]", "_", pathway.show)
  
  # ---------------------------------------------------------
  # TYPE 1: Hierarchy Plot
  # ---------------------------------------------------------
  hier_path <- file.path(output_base_dir, "Hierarchy", paste0("Hierarchy_", safe_pathway_name, ".png"))
  png(filename = hier_path, width = 1200, height = 700, res = 130)
  par(mar = c(3, 3, 5, 3), xpd = TRUE)
  
  netVisual_aggregate(cellchat, signaling = pathway.show, vertex.receiver = vertex.receiver)
  title(main = paste(pathway.show, "Signaling Network (Hierarchy)"), cex.main = 1.5, line = 2)
  
  dev.off()
  
  # ---------------------------------------------------------
  # TYPE 2: Circle Plot
  # ---------------------------------------------------------
  circle_path <- file.path(output_base_dir, "Circle", paste0("Circle_", safe_pathway_name, ".png"))
  png(filename = circle_path, width = 900, height = 950, res = 130)
  par(mar = c(3, 3, 6, 3), xpd = TRUE)
  
  netVisual_aggregate(cellchat, signaling = pathway.show, layout = "circle",
                      vertex.label.cex = 0.8, margin = 0.2, 
                      arrow.size = 1.5, arrow.width = 1.5)
  title(main = paste(pathway.show, "Signaling Network (Circle)"), cex.main = 1.5, line = 3)
  
  dev.off()
  
  # ---------------------------------------------------------
  # TYPE 3: Chord Diagram
  # ---------------------------------------------------------
  chord_path <- file.path(output_base_dir, "Chord", paste0("Chord_", safe_pathway_name, ".png"))
  png(filename = chord_path, width = 900, height = 950, res = 130)
  par(mar = c(3, 3, 6, 3), xpd = TRUE)
  
  netVisual_aggregate(cellchat, signaling = pathway.show, layout = "chord")
  title(main = paste(pathway.show, "Signaling Network (Chord)"), cex.main = 1.5, line = 3)
  
  dev.off()
  
  # ---------------------------------------------------------
  # TYPE 4: Heatmap (Fixed with explicit print call)
  # ---------------------------------------------------------
  heatmap_path <- file.path(output_base_dir, "Heatmap", paste0("Heatmap_", safe_pathway_name, ".png"))
  
  tryCatch({
    png(filename = heatmap_path, width = 800, height = 800, res = 130)
    
    # FIX: Wrap the heatmap in print() so it forces R to draw it inside the loop
    p_heat <- netVisual_heatmap(cellchat, signaling = pathway.show, color.heatmap = "Reds")
    print(p_heat)
    
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off() 
    message(paste0("⚠️ Skipped Heatmap for '", pathway.show, "' - Reason: No distinct signaling variance to plot."))
  })
  
}

#-----------------------------------------------------------------------------------------------------------------------------------
library(dplyr)

# ==============================================================================
# 1. SETUP PATHS & VARIABLES
# ==============================================================================
output_base_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE/Signaling_Pathways"
all_pathways    <- cellchat@netP$pathways

# Define your target receiver cell indices for the hierarchy plot (e.g., clusters 1 to 4)
vertex.receiver <- seq(1, 4) 

# ==============================================================================
# 2. DYNAMIC CELL TYPE GROUPING (Tailored to your exact data)
# ==============================================================================
cluster_names <- levels(cellchat@idents)

group_mapping <- case_when(
  cluster_names == "Células Cancerígenas" ~ "Tumor",
  cluster_names == "Células Epiteliais"   ~ "Epithelial",
  cluster_names == "Células Endoteliais"  ~ "Endothelial",
  cluster_names == "Células Mielóides"    ~ "Myeloid",
  TRUE                                    ~ "Other"
)

group.cellType <- group_mapping
names(group.cellType) <- cluster_names

# ==============================================================================
# 3. CREATE DIRECTORY STRUCTURE
# ==============================================================================
modes <- c("Hierarchy", "Circle", "Chord", "Heatmap", "Chord_Cell_Grouped", 
           "Contribution", "Individual_Hierarchy", "Individual_Circle", "Individual_Chord")

for (mode in modes) {
  dir.create(file.path(output_base_dir, mode), recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 4. MASTER AUTOMATION LOOP
# ==============================================================================
for (pathway.show in all_pathways) {
  
  # Create a Windows-safe file name for the current pathway
  safe_pathway_name <- gsub("[^[:alnum:]_]", "_", pathway.show)
  
  message(paste("🔄 Processing pathway:", pathway.show))
  
  # ----------------------------------------------------------------------------
  # TYPE 1: Hierarchy Plot (Aggregated Pathway Level)
  # ----------------------------------------------------------------------------
  hier_path <- file.path(output_base_dir, "Hierarchy", paste0("Hierarchy_", safe_pathway_name, ".png"))
  png(filename = hier_path, width = 1200, height = 700, res = 130)
  par(mar = c(3, 3, 5, 3), xpd = TRUE)
  
  netVisual_aggregate(cellchat, signaling = pathway.show, vertex.receiver = vertex.receiver)
  title(main = paste(pathway.show, "Signaling Network (Hierarchy)"), cex.main = 1.5, line = 2)
  
  dev.off()
  
  # ----------------------------------------------------------------------------
  # TYPE 2: Circle Plot (Aggregated Pathway Level)
  # ----------------------------------------------------------------------------
  circle_path <- file.path(output_base_dir, "Circle", paste0("Circle_", safe_pathway_name, ".png"))
  png(filename = circle_path, width = 900, height = 950, res = 130)
  par(mar = c(3, 3, 6, 3), xpd = TRUE)
  
  netVisual_aggregate(cellchat, signaling = pathway.show, layout = "circle",
                      vertex.label.cex = 0.8, margin = 0.2, 
                      arrow.size = 1.5, arrow.width = 1.5)
  title(main = paste(pathway.show, "Signaling Network (Circle)"), cex.main = 1.5, line = 3)
  
  dev.off()
  
  # ----------------------------------------------------------------------------
  # TYPE 3: Chord Diagram (Aggregated Pathway Level)
  # ----------------------------------------------------------------------------
  chord_path <- file.path(output_base_dir, "Chord", paste0("Chord_", safe_pathway_name, ".png"))
  png(filename = chord_path, width = 900, height = 950, res = 130)
  par(mar = c(3, 3, 6, 3), xpd = TRUE)
  
  netVisual_aggregate(cellchat, signaling = pathway.show, layout = "chord")
  title(main = paste(pathway.show, "Signaling Network (Chord)"), cex.main = 1.5, line = 3)
  
  dev.off()
  
  # ----------------------------------------------------------------------------
  # TYPE 4: Heatmap (Aggregated Pathway Level - Fixed with explicit print)
  # ----------------------------------------------------------------------------
  heatmap_path <- file.path(output_base_dir, "Heatmap", paste0("Heatmap_", safe_pathway_name, ".png"))
  
  tryCatch({
    png(filename = heatmap_path, width = 800, height = 800, res = 130)
    p_heat <- netVisual_heatmap(cellchat, signaling = pathway.show, color.heatmap = "Reds")
    print(p_heat) 
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off() 
    message(paste0("  ⚠️ Skipped Heatmap for '", pathway.show, "' - Reason: No distinct signaling variance to plot."))
  })
  
  # ----------------------------------------------------------------------------
  # TYPE 5: Multi-Group Chord Diagram (Dynamic Cell Type Level)
  # ----------------------------------------------------------------------------
  chord_cell_path <- file.path(output_base_dir, "Chord_Cell_Grouped", paste0("Grouped_Chord_", safe_pathway_name, ".png"))
  
  tryCatch({
    png(filename = chord_cell_path, width = 1000, height = 1000, res = 130)
    netVisual_chord_cell(cellchat, signaling = pathway.show, group = group.cellType, 
                         title.name = paste0(pathway.show, " Grouped Signaling Network"))
    dev.off()
  }, error = function(e) { 
    if (dev.cur() > 1) dev.off()
    message(paste0("  ⚠️ Skipped Grouped Chord for '", pathway.show, "'")) 
  })
  
  # ----------------------------------------------------------------------------
  # TYPE 6: Ligand-Receptor Pair Contribution Bar Plot
  # ----------------------------------------------------------------------------
  contrib_path <- file.path(output_base_dir, "Contribution", paste0("Contribution_", safe_pathway_name, ".png"))
  png(filename = contrib_path, width = 800, height = 600, res = 130)
  
  p_contrib <- netAnalysis_contribution(cellchat, signaling = pathway.show)
  print(p_contrib) 
  
  dev.off()
  
  # ----------------------------------------------------------------------------
  # TYPE 7: INDIVIDUAL LIGAND-RECEPTOR PAIR PLOTS
  # ----------------------------------------------------------------------------
  pairLR.all <- extractEnrichedLR(cellchat, signaling = pathway.show, geneLR.return = FALSE)
  
  if (!is.null(pairLR.all) && nrow(pairLR.all) > 0) {
    
    # Grab the top-ranked L-R pair for individual visualization
    LR.show <- pairLR.all[1, ] 
    safe_LR_name <- gsub("[^[:alnum:]_]", "_", LR.show)
    
    # A. Individual Hierarchy
    ind_hier_path <- file.path(output_base_dir, "Individual_Hierarchy", paste0("Ind_Hierarchy_", safe_pathway_name, "_", safe_LR_name, ".png"))
    png(filename = ind_hier_path, width = 1200, height = 700, res = 130)
    par(mar = c(3, 3, 5, 3), xpd = TRUE)
    netVisual_individual(cellchat, signaling = pathway.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver)
    title(main = paste("L-R Pair:", LR.show, "(Hierarchy)"), cex.main = 1.2, line = 2)
    dev.off()
    
    # B. Individual Circle
    ind_circle_path <- file.path(output_base_dir, "Individual_Circle", paste0("Ind_Circle_", safe_pathway_name, "_", safe_LR_name, ".png"))
    png(filename = ind_circle_path, width = 900, height = 950, res = 130)
    par(mar = c(3, 3, 6, 3), xpd = TRUE)
    netVisual_individual(cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle",
                         vertex.label.cex = 0.8, margin = 0.2)
    title(main = paste("L-R Pair:", LR.show, "(Circle)"), cex.main = 1.2, line = 3)
    dev.off()
    
    # C. Individual Chord
    ind_chord_path <- file.path(output_base_dir, "Individual_Chord", paste0("Ind_Chord_", safe_pathway_name, "_", safe_LR_name, ".png"))
    
    tryCatch({
      png(filename = ind_chord_path, width = 900, height = 950, res = 130)
      netVisual_individual(cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")
      dev.off()
    }, error = function(e) { 
      if (dev.cur() > 1) dev.off()
      message(paste0("  ⚠️ Skipped Individual Chord for '", LR.show, "'")) 
    })
    
  } else {
    message(paste0("  ℹ️ No specific L-R pairs found for pathway: ", pathway.show))
  }
}

message("🎉 Done! All plots have been successfully generated and saved.")


#-----------------------------------------------------------------------------------------------------------------------------------





# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


#-----------------------------------------------------------------------------------------------------------------------------------


# 1. Define paths and make sure the new "Bubble" directory exists
output_base_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE/Signaling_Pathways"
dir.create(file.path(output_base_dir, "Bubble_Plots"), recursive = TRUE, showWarnings = FALSE)

# 2. Automatically get all 4 of your cell groups dynamically (indices 1 to 4)
all_clusters <- 1:length(levels(cellchat@idents)) 

# 3. Get all available pathways
all_pathways <- cellchat@netP$pathways

# ==============================================================================
# AUTOMATION LOOP FOR BUBBLE PLOTS
# ==============================================================================
for (pathway.show in all_pathways) {
  
  safe_pathway_name <- gsub("[^[:alnum:]_]", "_", pathway.show)
  bubble_path <- file.path(output_base_dir, "Bubble_Plots", paste0("Bubble_", safe_pathway_name, ".png"))
  
  # netVisual_bubble outputs a ggplot object. If it fails due to no matching pairs, 
  # tryCatch keeps the script running smoothly.
  tryCatch({
    # We increase the canvas size slightly because bubble plots can have long text labels
    png(filename = bubble_path, width = 1100, height = 900, res = 130)
    
    # Generate the bubble plot across ALL your cell groups for this specific pathway
    p_bubble <- netVisual_bubble(cellchat, 
                                 sources.use = all_clusters, 
                                 targets.use = all_clusters, 
                                 signaling = pathway.show, 
                                 remove.isolate = TRUE) # TRUE hides uninformative, blank rows
    
    # Force R to print the ggplot object to the png device
    print(p_bubble)
    
    dev.off()
    
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message(paste0("  ⚠️ Skipped Bubble plot for '", pathway.show, "' - Reason: No active L-R pairs found for this layout."))
  })
}
#-----------------------------------------------------------------------------------------------------------------------------------



#ABAIXO SÃO EXEMPLOS DO SITE, ADAPTAR PARA O CASO QUE ESTIVER SENDO AVALIADO

# set the order of interacting cell pairs on x-axis
# (4) Default: first sort cell pairs based on the appearance of sources in levels(object@idents), and then based on the appearance of targets in levels(object@idents)
# (5) sort cell pairs based on the targets.use defined by users
netVisual_bubble(cellchat, targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.target = T)
# (6) sort cell pairs based on the sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T)
# (7) sort cell pairs based on the sources.use and then targets.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T)
# (8) sort cell pairs based on the targets.use and then sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)

#(B) Chord diagram
#Similar to Bubble plot, CellChat provides a function netVisual_chord_gene for drawing Chord diagram to

#show all the interactions (L-R pairs or signaling pathways) from some cell groups to other cell groups. Two special cases: one is showing all the interactions sending from one cell groups and the other is showing all the interactions received by one cell group.

#show the interactions inputted by USERS or certain signaling pathways defined by USERS

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 6, targets.use = c(1:3), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3), targets.use = 6, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3), targets.use = c(1:8), signaling = c("CCL","CXCL"),legend.pos.x = 8)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(1:8), slot.name = "netP", legend.pos.x = 10)



#-----------------------------------------------------------------------------------------------------------------------------------



# 1. Define paths and create the directory
output_base_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE/Signaling_Pathways"
dir.create(file.path(output_base_dir, "Chord_Gene_Level"), recursive = TRUE, showWarnings = FALSE)

# 2. Automatically get all 4 of your cell groups dynamically
all_clusters <- 1:length(levels(cellchat@idents)) 

# 3. Get all available pathways
all_pathways <- cellchat@netP$pathways

# ==============================================================================
# AUTOMATION LOOP FOR GENE-LEVEL CHORD DIAGRAMS (With Larger Text)
# ==============================================================================
for (pathway.show in all_pathways) {
  
  safe_pathway_name <- gsub("[^[:alnum:]_]", "_", pathway.show)
  chord_gene_path <- file.path(output_base_dir, "Chord_Gene_Level", paste0("Chord_Gene_", safe_pathway_name, ".png"))
  
  tryCatch({
    png(filename = chord_gene_path, width = 1400, height = 1400, res = 140)
    
    # We add a slight top margin cushion via base R par() to protect the larger title text
    par(mar = c(2, 2, 6, 2), xpd = TRUE)
    
    netVisual_chord_gene(cellchat, 
                         sources.use = all_clusters, 
                         targets.use = all_clusters, 
                         signaling = pathway.show,
                         lab.cex = 1.2,            # FIX 1: Increased gene text size (Up from 0.5)
                         small.gap = 1,            
                         big.gap = 8,              
                         legend.pos.y = 30)       
    
    # FIX 2: Manually overlay a crisp, large title at the top
    title(main = paste(pathway.show, "Gene-Level Signaling Network"), 
          cex.main = 2.0,   # Doubles the title text size
          line = 3)         # Pushes it up safely above the circle
    
    dev.off()
    
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message(paste0("  ❌ Truly skipped '", pathway.show, "' - Actual Error: ", e$message))
  })
}


#-----------------------------------------------------------------------------------------------------------------------------------



#IMPORTANTE, SE HOUVER OVERLAP:
#NB: Please ignore the note when generating the plot such as “Note: The first link end is drawn out of sector ‘MIF’.”.
#If the gene names are overlapped, you can adjust the argument small.gap by decreasing the value.

#Plot the signaling gene expression distribution using violin/dot plot

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")
print(as.numeric(execution.time, units = "secs"))

#By default, plotGeneExpression only shows the expression of signaling genes related to the inferred significant communications. USERS can show the expression of all signaling genes related to one signaling pathway by

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)
execution.time = Sys.time() - ptm

### PARTE 4: ANÁLISE DE SISTEMAS DA REDE DE COMUNICAÇÃO CÉLULA-CÉLULA

#>To facilitate the interpretation of the complex intercellular communication networks, CellChat quantitively measures networks through methods abstracted from graph theory, pattern recognition and manifold learning.
#>
#Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling

#A) Compute and visualize the network centrality scores
ptm = Sys.time()
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

########## LOOP QUE EU FIZ #####################################################
# Defina o diretório onde os arquivos serão salvos
output_dir <- "C:/Users/clayt/OneDrive/Documentos/CellCommunication/TNBC/NE/senders_and_receivers"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Medir tempo de início
ptm <- Sys.time()

# Calcular centralidade
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Loop por via de sinalização
for (i in 1:length(pathways.show.all)) {
  
  pathway_name <- pathways.show.all[i]
  message("Processando: ", pathway_name)
  
  ## -------------------------------
  ## 1. Heatmap de Centralidade
  ## -------------------------------
  pdf(
    file = file.path(output_dir, paste0(pathway_name, "_centrality_heatmap.pdf")),
    width = 8,
    height = 2.5
  )
  netAnalysis_signalingRole_network(
    cellchat,
    signaling = pathway_name,
    width = 8,
    height = 2.5,
    font.size = 10
  )
  dev.off()
  
  ## -------------------------------
  ## 2. Gráfico de Cordas (Chord Diagram)
  ## -------------------------------
  # Subset da via
  cellchat_sub <- subsetCommunication(cellchat, signaling = pathway_name, slot.name = "netP")
  
  # Ajustar altura do gráfico com base no número de interações
  n_links <- nrow(cellchat_sub)
  plot_height <- max(6, min(12, 0.25 * n_links))  # Altura entre 6 e 12, proporcional
  
  pdf(
    file = file.path(output_dir, paste0(pathway_name, "_chord_diagram.pdf")),
    width = 8,
    height = plot_height
  )
  netVisual_aggregate(
    cellchat,
    signaling = pathway_name,
    layout = "chord",
    slot.name = "netP"
  )
  dev.off()
}

# Tempo final
cat("Todos os gráficos foram gerados e salvos em:\n", output_dir, "\n")
cat("Tempo total de execução:\n")
print(Sys.time() - ptm)

################################################################################

#(B) Visualize dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#(C) Identify signals contributing the most to outgoing or incoming signaling of certain cell groups

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.new()
# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))
ht



#-----------------------------------------------------------------------------------------------------------------------------------



library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

# ==============================================================================
# 1. SETUP PATHS & VARIABLES
# ==============================================================================
output_base_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE/Signaling_Pathways"
all_pathways    <- cellchat@netP$pathways

# Create a dedicated directory for the Systems/Role Analysis outputs
role_dir <- file.path(output_base_dir, "Systems_Role_Analysis")
dir.create(role_dir, recursive = TRUE, showWarnings = FALSE)

message("🧮 Computing network centrality scores (this might take a moment)...")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# ==============================================================================
# PART A: GLOBAL 2D SIGNALING ROLES (Scatter Plots)
# ==============================================================================
message("📊 Generating 2D Signaling Role Scatter Plots...")
scatter_path <- file.path(role_dir, "Global_Signaling_Role_Scatter.png")

# Open device
png(filename = scatter_path, width = 1400, height = 700, res = 130)

# gg1 looks at all pathways aggregated; gg2 tracks your specific discovered pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat) + ggtitle("All Pathways Aggregated")
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = all_pathways) + ggtitle("Discovered Pathways Summary")

# Combine them side by side and print
print(gg1 + gg2)

dev.off()

# ==============================================================================
# PART B: GLOBAL OUTGOING VS INCOMING HEATMAP PATTERNS
# ==============================================================================
message("🔥 Generating Global Outgoing vs Incoming Pattern Heatmaps...")
pattern_path <- file.path(role_dir, "Global_Pattern_Heatmaps.png")

png(filename = pattern_path, width = 1400, height = 900, res = 130)

# Stripped custom titles so CellChat uses its native layout settings
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

# Draw them side-by-side using ComplexHeatmap's engine
draw(ht1 + ht2)

dev.off()

# ==============================================================================
# PART C: PATHWAY-SPECIFIC ROLE HEATMAP
# ==============================================================================
message("🔥 Generating Pathway-Specific Role Heatmap Summary...")
spec_heatmap_path <- file.path(role_dir, "Pathway_Specific_Signaling_Role_Heatmap.png")

png(filename = spec_heatmap_path, width = 900, height = 900, res = 130)

# Stripped custom title here as well
ht_spec <- netAnalysis_signalingRole_heatmap(cellchat, signaling = all_pathways)
draw(ht_spec)

dev.off()

# ==============================================================================
# PART D: GENE EXPRESSION DISTRIBUTION (Violin Plots per Pathway)
# ==============================================================================
violin_dir <- file.path(output_base_dir, "Gene_Expression_Violins")
dir.create(violin_dir, recursive = TRUE, showWarnings = FALSE)

message("🎻 Generating Gene Expression Violin Plots...")
for (pathway.show in all_pathways) {
  
  safe_pathway_name <- gsub("[^[:alnum:]_]", "_", pathway.show)
  violin_path <- file.path(violin_dir, paste0("Expression_", safe_pathway_name, ".png"))
  
  tryCatch({
    png(filename = violin_path, width = 1000, height = 700, res = 130)
    
    p_vln <- plotGeneExpression(cellchat, signaling = pathway.show, enriched.only = TRUE, type = "violin")
    print(p_vln)
    
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message(paste0("  ⚠️ Skipped Violin Plot for '", pathway.show, "'"))
  })
}

message("🎉 Execution finished! All files are saved successfully in your OneDrive folder.")



#-----------------------------------------------------------------------------------------------------------------------------------



### Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
#In addition to exploring detailed communications for individual pathways, an important question is how multiple cell groups and signaling pathways coordinate to function. CellChat employs a pattern recognition method to identify the global communication patterns.

#(A) Identify and visualize outgoing communication pattern of secreting cells

#Outgoing patterns reveal how the sender cells (i.e. cells as signal source) coordinate with each other as well as how they coordinate with certain signaling pathways to drive communication.

library(NMF)
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 2/2
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare

library(ggalluvial)
#Here we run selectK to infer the number of patterns.
selectK(cellchat, pattern = "outgoing")

#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 6.
#verificar para o dataset quando tem o drop. Aqui foi 4
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.new()
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#(B) Identify and visualize incoming communication pattern of target cells

#Incoming patterns show how the target cells (i.e. cells as signal receivers) coordinate with each other as well as how they coordinate with certain signaling pathways to respond to incoming signals.
selectK(cellchat, pattern = "incoming")

#Cophenetic values begin to drop when the number of incoming patterns is 3.
#Veja qual é o drop number pra esse dataset.
dev.new()
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.new()
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

#### Manifold and classification learning analysis of signaling networks
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 147.8175

#SALVAR: adaptar local para salvar
saveRDS(cellchat, file = "C:/Users/clayt/OneDrive/Documentos/CellCommunication/TNBC/NE/cellchat_BRCA_Pembro_On_TNBC_NE.rds")
#terminar aqui conforme a vinheta:
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#part-iv-systems-analysis-of-cell-cell-communication-network



#-----------------------------------------------------------------------------------------------------------------------------------


library(NMF)
library(ggalluvial)
library(dplyr)

# ==============================================================================
# 1. PATHS & INITIAL PARAMETERS
# ==============================================================================
output_base_dir <- "D:/Onedrive/Renata_Proj/Test_Proj/check/RT_NE/Signaling_Pathways"
global_dir      <- file.path(output_base_dir, "Global_Coordination_Analysis")
dir.create(global_dir, recursive = TRUE, showWarnings = FALSE)

# Start execution timer
ptm <- Sys.time()

# ==============================================================================
# 2. AUTOMATED NMF PATTERN COUNT DETERMINATION (Fixed ggplot extraction)
# ==============================================================================
message("🤖 Calculating optimal NMF factorization parameters...")
k_range <- 2:7

# Run CellChat's mathematical rank estimations
out_k_plot <- selectK(cellchat, pattern = "outgoing", k.range = k_range)
in_k_plot  <- selectK(cellchat, pattern = "incoming", k.range = k_range)

# FIX: Extract the metric data frames directly from the ggplot objects
out_metrics <- out_k_plot$data %>% filter(Measure == "Cophenetic")
in_metrics  <- in_k_plot$data %>% filter(Measure == "Cophenetic")

# Automatically pick the best pattern counts by maximizing Cophenetic stability
best_out_k <- out_metrics$k[which.max(out_metrics$score)]
best_in_k  <- in_metrics$k[which.max(in_metrics$score)]

message(paste("👉 Automated Selection: Outgoing Patterns (k) =", best_out_k))
message(paste("👉 Automated Selection: Incoming Patterns (k) =", best_in_k))

# ==============================================================================
# 3. IDENTIFY & VISUALIZE OUTGOING COMMUNICATION PATTERNS
# ==============================================================================
message("🌊 Processing Outgoing Patterns...")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = best_out_k)

# Save Outgoing River Plot
png(filename = file.path(global_dir, "Patterns_Outgoing_River.png"), width = 1100, height = 900, res = 130)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

# Save Outgoing Dot Plot 
png(filename = file.path(global_dir, "Patterns_Outgoing_Dot.png"), width = 1000, height = 800, res = 130)
print(netAnalysis_dot(cellchat, pattern = "outgoing"))
dev.off()

# ==============================================================================
# 4. IDENTIFY & VISUALIZE INCOMING COMMUNICATION PATTERNS
# ==============================================================================
message("🌊 Processing Incoming Patterns...")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = best_in_k)

# Save Incoming River Plot
png(filename = file.path(global_dir, "Patterns_Incoming_River.png"), width = 1100, height = 900, res = 130)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

# Save Incoming Dot Plot
png(filename = file.path(global_dir, "Patterns_Incoming_Dot.png"), width = 1000, height = 800, res = 130)
print(netAnalysis_dot(cellchat, pattern = "incoming"))
dev.off()

# ==============================================================================
# 5. MANIFOLD EMBEDDINGS (Functional vs Structural Similarity)
# ==============================================================================
message("🧬 Clustering Networks by Functional Similarity Space...")
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Save Functional Embedding
png(filename = file.path(global_dir, "Embedding_Functional_Global.png"), width = 900, height = 900, res = 130)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

png(filename = file.path(global_dir, "Embedding_Functional_ZoomIn.png"), width = 1200, height = 600, res = 130)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()

message("🧬 Clustering Networks by Structural Topology Similarity Space...")
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

# Save Structural Embedding
png(filename = file.path(global_dir, "Embedding_Structural_Global.png"), width = 900, height = 900, res = 130)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()

png(filename = file.path(global_dir, "Embedding_Structural_ZoomIn.png"), width = 1200, height = 600, res = 130)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

# ==============================================================================
# 6. SAVE COMPLETED CELLCHAT ANALYSIS OBJECT
# ==============================================================================
message("💾 Archiving updated CellChat data asset...")
saveRDS(cellchat, file = file.path(output_base_dir, "cellchat_completed_analysis.rds"))

# Final time tracking report
execution.time <- Sys.time() - ptm
cat("\n==================================================\n")
cat("🎉 Task Complete! All pattern graphics have been written.\n")
cat("Total execution time:\n")
print(execution.time)
cat("==================================================\n")