
seurat_obj=readRDS("/Users/yevhenakimov/Downloads/GSE261983_integrated2.rds")

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




# 1. Summarize per-sample metrics ----------------------------------
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
            cluster_label = paste0("C", cluster_rank))  # e.g. "C1", "C2", ...

# 2. Define the new x-axis order by ascending sd_fraction ----------
sample_order <- summary_df %>%
  arrange(sd_fraction) %>%
  pull(sample)

# 3. Re-factor the summary AND the full alluvial data --------------
# 3a) Summary table
summary_df <- summary_df %>%
  mutate(sample = factor(sample, levels = sample_order))

# 3b) Full data, with zero‐padding so every (sample, cluster) exists
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

# 4. Compute a scale‐factor for the dual‐axis in the top panel
sf <- 0.8*max(summary_df$total_cells) / max(summary_df$total_clusters)

# 5. Top panel: total cells (bars), clusters (line+points), SD (errorbar)
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
  # geom_errorbar(aes(ymin = 0,
  #                   ymax = sd_fraction * sf),
  #               color = "#e31a1c",
  #               width = 0.2) +
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

# 6. Bottom panel: the alluvial plot
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

# 7. Combine them, giving the top 25% of the height to the summary
final_plot <- top_plot / alluvial_plot + 
  plot_layout(heights = c(0.25, 0.8))

# final_plot <- top_plot+coord_flip() + alluvial_plot + 
#   plot_layout(heights = c(0.25, 0.8))

# 8. Draw
print(final_plot)














wide_incidence <- meta_full %>%
  dplyr::select(sample, bbknn_leiden,n) %>%
  tidyr::pivot_wider(
    names_from = bbknn_leiden,
    values_from = n,
    values_fill = 0)%>% column_to_rownames("sample") 

wide_incidence=as.data.frame(wide_incidence)
wide_incidence[wide_incidence!=0]=1

wide_incidence=wide_incidence[as.character(unique(meta_full$sample)),]

rownames(t(wide_incidence))
p=plot_circle_heatmap4(data=t(wide_incidence),shape_values=22,color = Greys5,)+theme_bw()+coord_fixed()






ggsave(paste0("/Users/yevhenakimov/Downloads/","wide_incidence.pdf"), p, width = 12, height = 7.5,device="pdf", dpi = 900)

dimnames(wide_incidence)


















