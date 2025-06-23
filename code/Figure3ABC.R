rm(list = ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(scales)
library(patchwork)
library(plyr)
library(Seurat)
library(ggthemes)

seurat_filename <- "/Users/yevhenakimov/Downloads/GSE261983_integrated2.rds"


local_save_base_umap             <- "base_umap.png"
local_save_sample_umap           <- "sample_umap.png"
local_save_cl28                  <-  "cl28.png"
local_save_cl26                  <- "cl26.png"
local_save_cl29                  <- "cl29.png"
local_save_pheno_input_contours  <-  "pheno_input_contours.png"
local_save_alluvium_stratum_plot <-  "alluvium_stratum_plot_c5.pdf"


# : Read data and source helper functions
seurat_obj <- readRDS(seurat_filename)
source("https://raw.githubusercontent.com/YevhenAkimov/graphics-R/main/graphics_functions.R")
source("https://raw.githubusercontent.com/YevhenAkimov/general_purpose_R/main/general_helpers.R")

# : Define color schemes
colorscheme_Miller_Stone    <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]][["Miller Stone"]]$value
Greys2                      <- c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525')
colorscheme_prism_6colors   <- c('#5F4690','#1D6996','#EDAD08','#CC503E','#94346E','#0F8554')

# : Base UMAP plot
base_umap <- ggscatter_colored(
  seurat_obj@reductions$umap@cell.embeddings,
  as.factor(seurat_obj@meta.data$bbknn_leiden),
  size_mult = 0.41,
  colors    = interpolate_colors(colorscheme_Miller_Stone, 30)
) +
  theme_void() +
  coord_fixed()

ggsave(
  filename = local_save_base_umap,
  plot     = base_umap,
  width    = 19, height = 19,
  dpi      = 600, device = "png"
)

# : Sample-colored UMAP plot
sample_umap <- ggscatter_colored(
  seurat_obj@reductions$umap@cell.embeddings,
  as.factor(seurat_obj@meta.data$sample),
  size_mult = 0.11,
  colors    = interpolate_colors(colorscheme_Miller_Stone, 30)
) +
  theme_void() +
  coord_fixed()

# : Single-cluster UMAP plots
single_cluster_plots <- list()
for (cl in unique(seurat_obj@meta.data$bbknn_leiden)) {
  cells_in_cl <- rownames(seurat_obj[[]])[seurat_obj@meta.data$bbknn_leiden == cl]
  single_cluster_plots[[as.character(cl)]] <- ggscatter_colored(
    seurat_obj@reductions$umap@cell.embeddings[cells_in_cl, ],
    as.factor(seurat_obj@meta.data[cells_in_cl, ]$sample),
    size_mult = 0.5,
    colors    = interpolate_colors(Greys2, 16)
  ) +
    theme_void() +
    coord_fixed()
}

# : Save selected single-cluster plots
ggsave(filename = local_save_cl28, plot = single_cluster_plots[["28"]] + theme(legend.position = "none"),
       width = 1.7, height = 1.7, dpi = 1000, device = "png")
ggsave(filename = local_save_cl26, plot = single_cluster_plots[["26"]] + theme(legend.position = "none"),
       width = 1.7, height = 1.7, dpi = 1000, device = "png")
ggsave(filename = local_save_cl29, plot = single_cluster_plots[["29"]] + theme(legend.position = "none"),
       width = 2.7, height = 2.7, dpi = 1000, device = "png")

# : Annotate sample names in metadata
df <- seurat_obj[[]]
df$row_names <- rownames(df)
seurat_obj@meta.data$sample_names <- sub("#.*", "", df$row_names)

# : Define marker list and calculate module scores
marker_list <- list(
  Oligodendrocytes = c("MOG",  "MBP",  "CNP",  "PLP1", "MAG"),
  Microglia        = c("PTPRC","CSF1R","AIF1","CX3CR1","P2RY12"),
  OPC              = c("PDGFRA","CSPG4","OLIG1","OLIG2","SOX10"),
  VIP_Interneurons = c("VIP",  "CALB2","HTR3A","RELN","TAC3")
)

seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = marker_list,
  name     = names(marker_list),
  ctrl     = 100,
  assay    = "RNA",
  verbose  = FALSE
)

# : Predict cell type assignment
seurat_obj$predicted_celltype <- apply(
  seurat_obj@meta.data[, c("VIP_Interneurons4", "OPC3", "Microglia2", "Oligodendrocytes1")],
  1,
  function(x) names(which.max(x))
)

# initialize all phenotype columns
for (ph in c("VIP_Interneurons", "OPC", "Microglia", "Oligodendrocytes")) {
  seurat_obj@meta.data[[ph]] <- 0
}
# assign scores to predicted types
seurat_obj@meta.data$VIP_Interneurons   [seurat_obj$predicted_celltype == "VIP_Interneurons4"]   <- seurat_obj$VIP_Interneurons4[seurat_obj$predicted_celltype == "VIP_Interneurons4"]
seurat_obj@meta.data$OPC                [seurat_obj$predicted_celltype == "OPC3"]                <- seurat_obj$OPC3[seurat_obj$predicted_celltype == "OPC3"]
seurat_obj@meta.data$Microglia          [seurat_obj$predicted_celltype == "Microglia2"]          <- seurat_obj$Microglia2[seurat_obj$predicted_celltype == "Microglia2"]
seurat_obj@meta.data$Oligodendrocytes  [seurat_obj$predicted_celltype == "Oligodendrocytes1"]  <- seurat_obj$Oligodendrocytes1[seurat_obj$predicted_celltype == "Oligodendrocytes1"]

# : Prepare inputs for phenotype contours
# Subsection: Prepare inputs for phenotype contours
coords_inp  <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
Pheno_input <- seurat_obj@meta.data[,c("VIP_Interneurons","OPC","Microglia","Oligodendrocytes")]

# create density list
dens_list <- lapply(seq_len(ncol(Pheno_input)), function(i){
  dens       <- get_density2d_spatstat(coords_inp, weights = Pheno_input[[i]], adjust = 0.1)
  dens$v[dens$v < 0] <- 0
  df         <- reshape2::melt(dens$v)
  df$y       <- dens$yrow[df$Var1]
  df$x       <- dens$xcol[df$Var2]
  df$value   <- min_max_normalization(df$value, 0, 1)
  df
})

# Subsection: Generate phenotype contour plots (fixed error)
contour_breaks <- seq(0.25, 1, length.out = 5)

p_base <- ggplot() +
  geom_point(data = coords_inp, aes(x = umap_1, y = umap_2), size = 0.2, color = "#a1a1a1") +
  geom_point(data = coords_inp, aes(x = umap_1, y = umap_2), size = 0.11, color = "#f1f1f1") +
  theme_void() + coord_fixed()

# iterate over indices to avoid equality comparison error
p_final <- Reduce(function(p, idx){
  df_i     <- dens_list[[idx]]
  base_col= colorspace::darken(  colorscheme_prism_6colors[idx],0.2,  space="combined", method = "relative" )

  p +
    geom_contour_filled(data = df_i, aes(x = x, y = y, z = value, alpha = ..level..),
                        fill = base_col, breaks = contour_breaks) +
    geom_contour(data = df_i, aes(x = x, y = y, z = value),
                 color = "white", alpha = 0.4, size = 0.5,
                 breaks = contour_breaks) +
    scale_alpha_discrete(range = c(0.1, 0.35))
}, seq_along(dens_list), init = p_base)

ggsave(
  filename = local_save_pheno_input_contours,
  plot     = p_final + theme_void(),
  width    = 12, height = 12,
  dpi      = 1200,
  device   = "png"
)

# : Identify outliers via PCA and UMAP
smpl_mean_pca <- summarize_columns(as.data.frame(t(seurat_obj@assays$RNA$data)),
                                   seurat_obj@meta.data$sample_names,
                                   order_rows = FALSE, indices = FALSE)
prcomp_res    <- prcomp(smpl_mean_pca)
umap_res      <- uwot::umap(prcomp_res$x, n_neighbors = 4)
high=rownames(umap_res)
names(high)=rownames(umap_res)
ggscatter_colored(umap_res,umap_res,highlight_points =high)


# : Combined cluster highlighting
mult=0.6
plot_= add_ggpoint( as.data.frame(seurat_obj@reductions$umap@cell.embeddings), sizes=c(0.55)*mult,colors=c("#a1a1a1"),alphas=1 )+theme_void()+coord_fixed()
plot_= add_ggpoint( as.data.frame(seurat_obj@reductions$umap@cell.embeddings),plot_, sizes=c(0.40)*mult,colors=c("#e1e1e1"),alphas=1 )+theme_void()+coord_fixed()
plot_

p=list()
for (i in unique(seurat_obj@meta.data$bbknn_leiden)) {
  c_names=rownames(seurat_obj@meta.data[seurat_obj@meta.data$bbknn_leiden==i,])
  p[[as.character(i)]] =  add_ggpoint( seurat_obj@reductions$umap@cell.embeddings[c_names,], plot_ ,sizes=c(1.80)*mult,colors=c("#be2312"),alphas=1 )+coord_fixed()
  p[[as.character(i)]] =  add_ggpoint(  seurat_obj@reductions$umap@cell.embeddings[c_names,],p[[as.character(i)]] , sizes=c(0.50)*mult,colors=c("#f1f1f1"),alphas=1 )+coord_fixed()
}
p[["29"]]
p[["28"]]
p[["26"]]


# : Calculate sample–cluster fractions and alluvial plot

meta_=seurat_obj@meta.data

meta_$sample=as.factor(meta_$sample)
meta_$bbknn_leiden=as.factor(meta_$bbknn_leiden)

## find fraction of each sample  in each cluster
meta_fract_laiden <- meta_ %>%
  group_by(sample, bbknn_leiden) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(bbknn_leiden) %>%
  dplyr::mutate(fraction = n / sum(n)) %>%
  ungroup()
meta_fract <- meta_ %>%
  group_by(sample, bbknn_leiden) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  dplyr::mutate(fraction = n / sum(n)) %>%
  ungroup()




#  Summarize per-sample metrics 
summary_df <- meta_fract %>%
  group_by(sample) %>%
  dplyr::summarise(
    total_cells    = sum(n),
    total_clusters = n_distinct(bbknn_leiden),
    sd_fraction    = sd(fraction) 
  ) %>%
  ungroup()
cluster_sd <- meta_fract %>%
  group_by(bbknn_leiden) %>%
  dplyr::summarise(sd_cluster = sd(fraction)) %>%
  arrange(sd_cluster) %>%
  dplyr::mutate(cluster_rank = row_number())    # 1 = smallest SD

cluster_map <- cluster_sd %>%
  transmute(bbknn_leiden,
            cluster_label = paste0("C", cluster_rank)) 

# Define the new x-axis order by ascending sd_fraction 
sample_order <- summary_df %>%
  arrange(sd_fraction) %>%
  pull(sample)

#  Re-factor the summary AND the full alluvial data 
#  Summary table
summary_df <- summary_df %>%
  mutate(sample = factor(sample, levels = sample_order))

#  Full data, with zero‐padding so every (sample, cluster) exists
meta_full <- meta_fract %>%
  dplyr::mutate(
    sample         = factor(sample, levels = sample_order),
    bbknn_leiden   = factor(bbknn_leiden,
                            levels = sort(unique(bbknn_leiden)))
  ) %>%
  complete(
    sample, bbknn_leiden,
    fill = list(n = 0, fraction = 0)
  )

#  Compute a scale‐factor for the dual‐axis 
sf <- 0.8*max(summary_df$total_cells) / max(summary_df$total_clusters)

#  Top panel
top_plot <- ggplot(summary_df, aes(x = sample)) +
  geom_hline(yintercept = 0, color = "#313131",size=0.5)+
  geom_hline(yintercept =  30 * sf, color = "#661100",size=0.5, linetype="longdash")+
  geom_col(aes(y = total_cells), fill = "#c1c1c1",color="#717171",width=0.5,size=0.5,
  ) +
  geom_line(aes(y = total_clusters * sf, group = 1),
            color = "#1C1F33", size = 0.5) +
  geom_point(aes(y = total_clusters * sf),
             color = "#1C1F33", size = 2.5, shape=23, fill="#661100") +
  geom_hline(yintercept = 0, color = "#313131",size=0.5)+
  
  scale_y_continuous(
    name    = "Total cells",
    labels  = comma,
    sec.axis = sec_axis(~ . / sf,
                        name = "Clusters ")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    #axis.ticks.x  = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line()
  )


alluvial_plot <- ggplot(meta_full,
                        aes(x        = sample,
                            stratum  = bbknn_leiden,
                            alluvium = bbknn_leiden,
                            y        = fraction,
                            fill     = bbknn_leiden)) +
  geom_alluvium(stat          = "alluvium",
                lode.guidance = "forward",
                colour        = NA,
                alpha         = 0.5,
                width         = 0.2) +
  geom_stratum(width  = 0.5,
               colour = "#1C1F33",
               size   = 0.2) +
  scale_fill_manual(values = interpolate_colors((colorscheme_Miller_Stone),30)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Sample", y = "Percent of cells") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_line(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "right"
  )

final_plot <- top_plot / alluvial_plot + 
  plot_layout(heights = c(0.25, 0.8))
# : Save alluvial stratum plot
ggsave(
  filename = local_save_alluvium_stratum_plot,
  plot     = final_plot,
  width    = 6, height = 4,
  device   = "pdf",
  dpi      = 900
)