

rm(list=ls())
# ─────────────────────────────────────
library(arrow)     
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("https://raw.githubusercontent.com/YevhenAkimov/general_purpose_R/main/general_helpers.R")
source("https://raw.githubusercontent.com/YevhenAkimov/graphics-R/main/graphics_functions.R")

df <- read_parquet("/stability_long.parquet") ## output from prepare_cluster_stability_R.py


full_tbl   <- df %>% dplyr::filter(kind == "full")  
inter_tbl  <- df %>%  dplyr::filter(kind == "inter")  
assign_tbl <- df %>%  dplyr::filter(kind == "assign")


full_diag <- full_tbl %>%
  dplyr::filter(row == col) %>%
  dplyr::select(sample, iter, N, row, full_value = value)

matched_inter <- assign_tbl %>%         
  dplyr::select(sample, iter, N, row, col) %>% 
  left_join(inter_tbl,
            by = c("sample", "iter", "N", "row", "col")) %>%
  dplyr::rename(inter_value = value)

ratio_tbl <- matched_inter %>%
  left_join(full_diag,
            by = c("sample", "iter", "N", "row")) 


ratio_tbl_1=ratio_tbl
ratio_tbl_1=ratio_tbl_1[ratio_tbl_1$sample == "GSM8155485_RT00374N", ] 

ggplot()+
  geom_point(data = ratio_tbl_1, aes(x = full_value, y = inter_value)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Full vs Inter values",
       x = "Full value", y = "Inter value") +
  theme_bw()+facet_wrap(~ N, scales = "free")


zero_replace=(find_min_except_zero( ratio_tbl$full_value )+find_min_except_zero( ratio_tbl$inter_value ))/2


ratio_tbl$full_value = ratio_tbl$full_value  + zero_replace
ratio_tbl$full_value[ratio_tbl$full_value>1]=1

ratio_tbl$inter_value = ratio_tbl$inter_value + zero_replace
ratio_tbl$inter_value[ratio_tbl$inter_value>1]=1

ratio_tbl$ratio=ratio_tbl$inter_value/ ratio_tbl$full_value


ratio_tbl <- ratio_tbl %>%
  dplyr::mutate(ratio = inter_value / full_value) 



ratio_tbl$ratio[ratio_tbl$ratio >1] <- 1
ratio_tbl_col=reshape2::dcast(ratio_tbl,  N ~ sample, value.var = "col" , fun.aggregate=max)
ratio_tbl_row=reshape2::dcast(ratio_tbl,  N ~ sample, value.var = "row" , fun.aggregate=max)
ratio_tbl_clust_ratio=ratio_tbl_col/ratio_tbl_row
ratio_tbl_clust_ratio[ratio_tbl_clust_ratio>1]=1
ratio_tbl_mean=reshape2::dcast(ratio_tbl, N ~ sample, value.var = "ratio" , fun.aggregate=mean)
ratio_tbl_sd =reshape2::dcast(ratio_tbl , N ~ sample, value.var = "ratio" , fun.aggregate=sd, na.rm=T)


rownames(ratio_tbl_mean)=ratio_tbl_mean$N
rownames(ratio_tbl_clust_ratio)=ratio_tbl_mean$N
ratio_tbl_mean= ratio_tbl_mean[,-1]
ratio_tbl_clust_ratio= ratio_tbl_clust_ratio[,-1]

dim(ratio_tbl_clust_ratio)
ratio_tbl_mean=as.data.frame(ratio_tbl_mean)
ratio_tbl_clust_ratio=as.data.frame(ratio_tbl_clust_ratio)


substrs=c("390","387N","386N","383N","382N","378N","376N")

selected_samples <-colnames(ratio_tbl_mean)[sapply(colnames(ratio_tbl_mean), function(x) any(sapply(substrs, function(sub) grepl(sub, x))))]


ratio_tbl_clust_iter <- ratio_tbl %>%             
  dplyr::group_by(sample, N, iter) %>%
  dplyr::summarise(clust_ratio = max(col) / max(row), .groups = "drop")

ratio_tbl_clust_sd <- ratio_tbl_clust_iter %>%      
  dplyr::group_by(sample, N) %>%
  dplyr::summarise(sd_clust = sd(clust_ratio), .groups = "drop")

ratio_tbl_clust_sd_wide <- reshape2::dcast(        
  ratio_tbl_clust_sd, N ~ sample, value.var = "sd_clust"
)

rownames(ratio_tbl_clust_sd_wide) <- ratio_tbl_clust_sd_wide$N
ratio_tbl_clust_sd_wide          <- ratio_tbl_clust_sd_wide[, -1]

clust_sd_mean <- data.frame(
  N  = as.numeric(rownames(ratio_tbl_clust_sd_wide)),
  SD = rowMeans(ratio_tbl_clust_sd_wide, na.rm = TRUE)
) |>
  arrange(N) |>
  mutate(x_idx = row_number())

agree_sd_mean <- data.frame(
  N  = as.numeric(rownames(ratio_tbl_sd)),
  SD = rowMeans(ratio_tbl_sd[ , -1, drop = FALSE], na.rm = TRUE)  
) |>
  arrange(N) |>
  mutate(x_idx = row_number())

clust_sd_mean <- clust_sd_mean %>%
  arrange(N) %>%
  mutate(x_idx = row_number())

agree_sd_mean <- agree_sd_mean %>%
  arrange(N) %>%
  mutate(x_idx = row_number())




idx_breaks <- clust_sd_mean$x_idx
idx_labels <- clust_sd_mean$N
p_sd_clust <- ggplot(clust_sd_mean, aes(x = x_idx, y = SD)) +
  geom_area(fill = "#515151", alpha = 0.1) +       
  geom_line(color = "#121212", size = 1) +        
  geom_point(color = "#121212", size = 3,shape=18) +        
  scale_x_continuous(
    breaks = idx_breaks,
    labels = idx_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(title = "SD of detected species (%)", y = "SD") +
  theme_minimal() +
  theme(
    axis.ticks = element_line(),
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.x = element_blank()
  ) + geom_hline(yintercept = 0, linetype = "solid", color = "#a2a2a2",size=1)


p_sd_agree <- ggplot(agree_sd_mean, aes(x = x_idx, y = SD)) +
  geom_area(fill = "#515151", alpha = 0.1) +        
  geom_line(color = "#121212", size = 1) +     
  geom_point(color = "#121212", size =  3,shape=18) +    
  scale_x_continuous(
    breaks = idx_breaks,
    labels = idx_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(title = "SD of cluster agreement", y = "SD") +
  theme_minimal() +
  theme(
    axis.ticks = element_line(),
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.x = element_blank()
  )+ geom_hline(yintercept = 0, linetype = "solid", color = "#a2a2a2",size=1)

ratio_tbl_mean[ratio_tbl_mean<0.915]=0.915
ratio_tbl_clust_ratio[ratio_tbl_clust_ratio<0.8]=0.8


#### Heatmap plotting function for data
data         = t(ratio_tbl_clust_ratio[,selected_samples]) * 100
data_sizes   = t(ratio_tbl_mean[,selected_samples])
p_heat <- ggshape_heatmap(
  data   = data,
  data_sizes   = data_sizes,
  shape_values = 22,            
  size_range   = c(2, 7),
  colorscheme  = c('#9C0824','#BF1316','#D42922','#E96251','#EBA49A','#B0B0B0','#838383','#5D5D5D', '#3B3B3B' ,'#1E1E1E') ,
  value_label  = "Detected species (%)",
  size_label   = "Cluster agreement score",
  column_label  = "Cells sampled",
  row_label = "Sample ID",
  title        = "Cluster ratio heatmap",
  theme_choice = theme_minimal(),
  legend.key.size  = 0.6,
  legend.text.size = 10,
  na_rm        = FALSE,
  text.angle.x = 45,
  text.angle.y = 0,
  cluster_by   = "data",   # "data", "data_sizes", "both"
  cluster_rows = F,
  cluster_cols = F,
  darken_outline_color_coeff=-1,
  grid.pars    = list(
    grid.size     = 0.2,
    grid.color    = "grey80",
    grid.linetype = "dashed"
  )
)


final_plot <- p_sd_clust / p_sd_agree / p_heat +
  plot_layout(
    heights = c(1.3, 1.3, 6)
    #align   = "v"    # this makes the x-axes line up exactly
  ) &
  theme()

print(final_plot)


ggsave(paste0("/Users/yevhenakimov/out/ratio_heatmap/final_plot", Sys.Date(), ".pdf"),
       plot=final_plot, width =9, height = 6, dpi = 1300,device="pdf")




p=plot_circle_heatmap3( data_sizes=t(ratio_tbl_mean[,selected_samples]),data=t(ratio_tbl_clust_ratio[,selected_samples])*100,shape_values=22,colorscheme=Greys2.5,value_label="Detected species (%)",size_label="Cluster agrement",id1_label='Sample', id2_label="Cells sampled",title="")
p=p+theme_text1()+theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()
)+coord_fixed()

ggsave(paste0("/Users/yevhenakimov/out/ratio_heatmap/heatmap_ratio_", Sys.Date(), ".pdf"),
       plot=p, width =6, height = 5, dpi = 1300,device="pdf")
