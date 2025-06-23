input_rds_path <- ""
output_results_rds_path <- ""
output_plot_path <- ""

required_packages <- c("Seurat","dplyr","purrr","tidyr","ggplot2","tibble")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Install it with install.packages('%s').", pkg, pkg), call. = FALSE)
  }
}

library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

if (input_rds_path == "") {
  stop("Please set 'input_rds_path' before running.")
}

sce <- readRDS(input_rds_path)
meta <- sce@meta.data

calc_chao2 <- function(binary_mat){
  species_counts <- colSums(binary_mat)
  Q1 <- sum(species_counts == 1)
  Q2 <- sum(species_counts == 2)
  S_obs <- length(species_counts)
  if(Q2 > 0){
    S_chao2 <- S_obs + (Q1^2)/(2*Q2)
  } else {
    S_chao2 <- S_obs
  }
  c(S_chao2 = S_chao2, S_obs = S_obs, Q1 = Q1, Q2 = Q2)
}

set.seed(22)
all_samples <- unique(meta$sample)
n_boot <- 5000
nsamps <- seq(3, length(all_samples))

results <- map_dfr(nsamps, function(n_samp){
  message(sprintf("Processing sample size %2d / %2d …", n_samp, max(nsamps)))
  stats <- replicate(n_boot, {
    sel <- sample(all_samples, n_samp, replace = FALSE)
    submd <- filter(meta, sample %in% sel)
    ctab <- table(submd$sample, submd$bbknn_leiden)
    binmat <- (ctab > 0)*1L
    ch <- calc_chao2(binmat)
    c(obs = unname(ch["S_obs"]), chao2 = unname(ch["S_chao2"]))
  })
  if (is.null(dim(stats))) {
    stats <- matrix(stats, nrow = 2, byrow = TRUE, dimnames = list(c("obs","chao2"), NULL))
  }
  obs_vals <- stats[1,]
  chao2_vals <- stats[2,]
  tibble(sample_size = n_samp,
         obs_mean = mean(obs_vals),
         obs_sd = sd(obs_vals),
         chao2_mean = mean(chao2_vals),
         chao2_sd = sd(chao2_vals))
})

if (output_results_rds_path == "") {
  stop("Please set 'output_results_rds_path' before running.")
}

saveRDS(results, output_results_rds_path)

p <- ggplot(results, aes(x = sample_size)) +
  geom_ribbon(aes(ymin = chao2_mean - chao2_sd, ymax = chao2_mean + chao2_sd, fill = "Chao2 estimate"), alpha = 0.12, color = NA) +
  geom_line(aes(y = chao2_mean + chao2_sd), color = "#bd6a52", size = 0.3, show.legend = FALSE, alpha = 0.51) +
  geom_line(aes(y = chao2_mean - chao2_sd), color = "#bd6a52", size = 0.3, show.legend = FALSE, alpha = 0.51) +
  geom_ribbon(aes(ymin = obs_mean - obs_sd, ymax = obs_mean + obs_sd, fill = "Observed clusters"), alpha = 0.12, color = NA) +
  geom_line(aes(y = obs_mean + obs_sd), color = "#383838", size = 0.3, show.legend = FALSE, alpha = 0.51) +
  geom_line(aes(y = obs_mean - obs_sd), color = "#383838", size = 0.3, show.legend = FALSE, alpha = 0.51) +
  geom_line(aes(y = obs_mean, color = "Observed clusters"), size = 1.15, linetype = "solid") +
  geom_line(aes(y = chao2_mean, color = "Chao2 estimate"), size = 1.15, linetype = "solid") +
  scale_color_manual("", values = c("Observed clusters" = "#383838", "Chao2 estimate" = "#bd6a52")) +
  scale_fill_manual("", values = c("Observed clusters" = "#383838", "Chao2 estimate" = "#bd6a52")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(3, 17, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(22, 34), breaks = seq(22, 34, 2), labels = seq(22, 34, 2)) +
  labs(x = "Number of Samples", y = "Number of Clusters", title = "Stability of Chao2 Completeness Estimation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major = element_line(size = 0.2))

if (output_plot_path == "") {
  stop("Please set 'output_plot_path' before running.")
}

ggsave(output_plot_path, plot = p, width = 2.3, height = 4.2, device = "pdf")