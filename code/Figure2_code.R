library(dplyr)
library(MASS)
library(corrplot)
library(UpSetR)
library(ggplot2)
library(openxlsx)
library(knitr)
library(kableExtra)
library(gridExtra)
library(stats)
library(fitdistrplus)
library(openxlsx)
library(tidyr)
library(readxl)
library(reshape2)
library(fossil)
library(paletteer)
library(matrixcalc)
options(scipen=10)

#WD
setwd("/home/mitro/Downloads/Unique celltypes/Unique_celltypes/data/")

# Define file paths and organ names
file_paths <- list(
  Eye = "Eye - Tabula Sapiens.xlsx",
  Bladder = "Bladder - Tabula Sapiens.xlsx",
  Heart = "Heart - Tabula Sapiens.xlsx",
  Kidney = "Kidney - Tabula Sapiens.xlsx",
  Large_intestine = "Large intestine - Tabula Sapiens.xlsx",
  Bone_marrow = "Bone marrow - Tabula Sapiens.xlsx",
  Liver = "Liver - Tabula Sapiens.xlsx",
  Salivary_gland = "Salivary gland - Tabula Sapiens.xlsx",
  Lung = "Lung - Tabula Sapiens.xlsx",
  Lymph_node = "Lymph node - Tabula Sapiens.xlsx",
  Mammary_gland = "Mammary gland - Tabula Sapiens.xlsx",
  Pancreas = "Pancreas - Tabula Sapiens.xlsx",
  Prostate = "Prostate - Tabula Sapiens.xlsx",
  Skin = "Skin - Tabula Sapiens.xlsx",
  Small_intestine = "Small intestine - Tabula Sapiens.xlsx",
  Spleen = "Spleen - Tabula Sapiens.xlsx",
  Thymus = "Thymus - Tabula Sapiens.xlsx",
  Tongue = "Tongue - Tabula Sapiens.xlsx",
  Trachea = "Trachea - Tabula Sapiens.xlsx",
  Uterus = "Uterus - Tabula Sapiens.xlsx"
)

# Function to generate a positive definite covariance matrix
generate_cov_matrix <- function(num_organs, mean, sd) {
  # Step 1: Generate a random matrix with normal distribution
  random_matrix <- matrix(rnorm(num_organs^2, mean, sd), nrow = num_organs)
  
  # Step 2: Symmetrize the matrix to make it symmetric
  symmetric_matrix <- (random_matrix + t(random_matrix)) / 2
  
  # Step 3: Ensure positive definiteness by adding a constant to the diagonal
  diag(symmetric_matrix) <- diag(symmetric_matrix) + num_organs
  
  # Step 4: Check if the matrix is positive definite
  eigenvalues <- eigen(symmetric_matrix)$values
  if (all(eigenvalues > 0)) {
    print("Matrix is positive definite.")
  } else {
    print("Matrix is not positive definite. Adjusting the diagonal elements...")
    # Step 5: Make it positive definite by adding a small constant to the diagonal elements
    diag(symmetric_matrix) <- diag(symmetric_matrix) + abs(min(eigenvalues)) + 1e-6
  }
  
  return(symmetric_matrix)
}

# Initialize an empty data frame to store combined data
combined_data <- data.frame()

# Loop over the file paths and read the data
for (organ in names(file_paths)) {
  file_path <- file_paths[[organ]]
  # Read the data from the Excel file
  organ_data <- read_excel(file_path)
  # Add the "Organ" column
  organ_data <- organ_data %>% mutate(Organ = organ)
  # Append to the combined data frame
  combined_data <- bind_rows(combined_data, organ_data)
}

# Display the combined data
print(combined_data)

# Optionally, save the combined data to a new Excel file
write.csv(combined_data, "combined_cell_type_data_TS.csv", row.names = FALSE)

# Pivot the combined data to a wide format
data_wide <- combined_data %>%
  pivot_wider(names_from = Organ, values_from = Count, values_fill = list(Count = 0))

# Set row names as cell types
data_matrix <- as.data.frame(data_wide)
rownames(data_matrix) <- data_matrix$`Cell type`
data_matrix <- data_matrix %>% dplyr::select(-`Cell type`)

# Display the matrix
print(data_matrix)

# Calculate the correlation matrix
correlation_matrix <- cor(data_matrix, use = "pairwise.complete.obs")

corrplot::corrplot(correlation_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)
ggsave("data_correlation.svg", scale = 1.5, device = "svg")

#Simulate data
library(MASS)
library(corrplot)
library(ggplot2)

# Compute target covariance and correlation from the original data
target_cov_matrix <- cov(data_matrix)  
target_cor_matrix <- cor(data_matrix)

# Compute log-transformed data
epsilon <- 1e-5
data_matrix[data_matrix == 0] <- epsilon
log_data_matrix <- log(data_matrix)

# Compute means and covariance in log-space
log_means <- colMeans(log_data_matrix)
log_cov_matrix <- cov(log_data_matrix)
log_cor_matrix <- cor(log_data_matrix)

# Function to adjust covariance matrix in normal space
adjust_covariance_matrix <- function(target_cov, target_cor) {
  num_vars <- ncol(target_cov)
  adjusted_cov <- matrix(0, nrow = num_vars, ncol = num_vars)
  
  # Parameters of the best-fitting normal distribution for non-diagonal elements
  mean_norm <- 6.448077
  sd_norm <- 4.967736
  
  for (i in 1:num_vars) {
    for (j in 1:num_vars) {
      if (i == j) {
        # Diagonal elements (variances) remain unchanged except for transformation
        adjusted_cov[i, j] <- log(1 + max(0, target_cov[i, j] / (exp(log_cov_matrix[i, j]) - 1)))
      } else {
        # Sample from the fitted normal distribution
        sampled_cov <- rnorm(1, mean = mean_norm, sd = sd_norm)
        
        # Scale the sampled covariance based on correlation structure
        term <- target_cor[i, j] * sqrt(max(0, (exp(log_cov_matrix[i, i]) - 1) * (exp(log_cov_matrix[j, j]) - 1)))
        adjusted_cov[i, j] <- log(1 + max(0, sampled_cov * term))
      }
    }
  }
  
  # Ensure symmetry
  adjusted_cov <- (adjusted_cov + t(adjusted_cov)) / 2
  
  # Ensure positive definiteness
  library(Matrix)
  if (!is.positive.definite(adjusted_cov)) {
    adjusted_cov <- as.matrix(nearPD(adjusted_cov)$mat)
  }
  
  return(adjusted_cov)
}

# Adjust the normal-space covariance matrix
adjusted_log_cov_matrix <- adjust_covariance_matrix(target_cov_matrix, target_cor_matrix)

# Ensure no NaN or Inf values
adjusted_log_cov_matrix[is.na(adjusted_log_cov_matrix)] <- 0
adjusted_log_cov_matrix[is.infinite(adjusted_log_cov_matrix)] <- max(adjusted_log_cov_matrix[is.finite(adjusted_log_cov_matrix)])

# Ensure positive definiteness
library(Matrix)
if (!is.positive.definite(adjusted_log_cov_matrix)) {
  adjusted_log_cov_matrix <- as.matrix(nearPD(adjusted_log_cov_matrix)$mat)
}

# Increase small variances
#diag(adjusted_log_cov_matrix) <- diag(adjusted_log_cov_matrix) + 1e-6

#Generate log_means
sim_log_means <- rnorm(20, -9, 0.5)

# Generate multivariate normal data using the adjusted covariance
num_samples <- 1500
simulated_normal_data <- mvrnorm(n = num_samples, mu = sim_log_means, Sigma = adjusted_log_cov_matrix)

# Apply exponentiation to obtain lognormal data
simulated_lognormal_data <- exp(simulated_normal_data)

# Compute correlation matrix for simulated lognormal data
simulated_cor_matrix <- cor(simulated_lognormal_data)

# Compare the target vs. simulated correlation matrices
par(mfrow = c(1, 2))  # Set up side-by-side plots

corrplot(target_cor_matrix, method = "color", type = "upper",
         main = "Target Correlation (data_matrix)", tl.col = "black", mar = c(2, 2, 2, 2))

corrplot(simulated_cor_matrix, method = "color", type = "upper",
         main = "Simulated Lognormal Correlation", tl.col = "black", mar = c(2, 2, 2, 2))

# Reset plotting layout
par(mfrow = c(1, 1))

# Compute correlation difference for evaluation
cor_diff <- abs(target_cor_matrix - simulated_cor_matrix)

# Summary of correlation differences
summary(cor_diff)

#Plot simulated correlation
corrplot::corrplot(simulated_cor_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)

# Extract upper triangular values (excluding diagonal) for comparison
target_cor_values <- target_cor_matrix[upper.tri(target_cor_matrix , diag = FALSE)]
simulated_cor_values <- cor(simulated_lognormal_data)[upper.tri(target_cor_matrix, diag = FALSE)]

# Combine into a data frame
cor_comparison <- data.frame(
  Value = c(target_cor_values, simulated_cor_values),
  Type = rep(c("Tabula Sapiens Covariance", "Simulated Covariance"), each = length(target_cor_values))
)

ggplot(cor_comparison, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.5) +  # Transparency for overlapping curves
  labs(title = "Comparison of Target vs Simulated Covariance Distributions",
       x = "Pearson correlation", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Tabula Sapiens Covariance" = "blue", "Simulated Covariance" = "red"))

ggsave("sim_vs_target_correlation.svg", scale = 1.5, device = "svg")

#Other plots
# Convert the scaled data to integers and set very small values to zero
integer_simulated_data <- round(simulated_lognormal_data)
integer_simulated_data[integer_simulated_data < epsilon] <- 0

total_cells_in_20_organs <- 5000000000000

scaling_factor <- total_cells_in_20_organs / sum(colSums(integer_simulated_data))
scaled_dataframe <- integer_simulated_data * scaling_factor

scaled_dataframe <- scaled_dataframe[rowSums(scaled_dataframe != 0) > 0, ]


# Step 1: Calculate the row sums of your scaled dataframe
row_sums <- rowSums(scaled_dataframe)

# Step 2: Create a data frame for plotting
abundance_data <- data.frame(
  `Cell type` = as.factor(1:nrow(scaled_dataframe)),  # Row index as Cell type
  Count = row_sums  # Corresponding row sums
)

colnames(abundance_data) <- c("Cell type", "Count")

# Step 3: Plot using ggplot
ggplot(data = abundance_data, aes(x = reorder(`Cell type`, -Count), y = Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Amount of Cell Types from Tabula Sapiens", x = "Cell type", y = "Amount") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_blank(),
    legend
  )

abundance_data <- abundance_data[abundance_data$Count > 0, ]
nrow(abundance_data)

#Plot Fig2 C
# Calculate total cell types per organ
scaled_dataframe <- as.data.frame(scaled_dataframe)
colnames(scaled_dataframe) <- 1:20

scaled_dataframe_transformed <- scaled_dataframe %>%
  # Use mutate_all to apply the transformation to every column
  mutate_all(~ ifelse(. >= 0.5, 1, 0))

organ_cell_counts <- scaled_dataframe_transformed %>%
  summarise(across(everything(), ~ sum(. == 1)))  # Calculate total cell types for each organ

# Transpose the data and convert to data frame
organ_cell_counts <- as.data.frame(t(organ_cell_counts))
colnames(organ_cell_counts) <- "total"
organ_cell_counts$organ <- rownames(organ_cell_counts)

unique_counts <- sapply(scaled_dataframe_transformed, function(col) {
  sum(col == 1 & rowSums(scaled_dataframe_transformed) == 1)
})

# Assuming organ_cell_counts is already created and needs to be updated with 'unique' counts
organ_cell_counts$unique <- unique_counts

# Calculate shared cell types
organ_cell_counts$shared <- organ_cell_counts$total - organ_cell_counts$unique

# Convert to long format for plotting
plot_data <- organ_cell_counts %>%
  dplyr::select(organ, unique, shared) %>%
  pivot_longer(cols = c("unique", "shared"), names_to = "type", values_to = "count")

# Sort organs by total cell types
plot_data$organ <- reorder(plot_data$organ, -organ_cell_counts$total[match(plot_data$organ, organ_cell_counts$organ)])

ggplot(plot_data, aes(x = organ, y = count, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("unique" = "blue", "shared" = "lightblue")) +
  labs(
    title = "Total and Unique Cell Types Across Organs",
    x = "Organ",
    y = "Number of Cell Types",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),  # Black x-axis labels
    axis.text.y = element_text(size = 16, colour = "black"),               # Black y-axis labels
    text = element_text(size = 16)                                        # Overall text size
  ) +
  coord_flip()  # Flips the x and y axes for horizontal bars

ggsave("Sim_total_and_unique_types.svg", scale = 1.5, device = "svg")

#Estimation for 3-organs
oneorgan <- as.data.frame(scaled_dataframe$`9`)
colnames(oneorgan) <- "Count"
oneorgan <- oneorgan[oneorgan$Count != 0, ]
oneorgan <- as.data.frame(oneorgan)
colnames(oneorgan) <- "Count"
oneorgan <- round(oneorgan)
oneorgan$"Cell type" <- 1:nrow(oneorgan)
oneorgan$"prob" <- oneorgan$Count / sum(oneorgan$Count)
abundance_data_for_estimation <- oneorgan

n_samples <- 6
sample_fraction <- c(1e-11, 1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04)

#Abundance dataframe of the population
abundance_data_for_estimation$Count <- round(abundance_data_for_estimation$Count)
abundance_data_for_estimation$prob <- abundance_data_for_estimation$Count / sum(abundance_data_for_estimation$Count)

generate_freq_table <- function(sample_data, sample_num) {
  observed_freq <- table(sample_data)
  observed_freq_df <- as.data.frame(observed_freq)
  
  # Ensure columns are named correctly
  colnames(observed_freq_df) <- c("Value", paste0("Frequency_Sample_", sample_num))
  
  return(observed_freq_df)
}

# Define the column names
columns <- c("Sample coverage", "Found in sample", "ICE", "Chao2", "Jack1", "Jack2")

# Create an empty data frame with these columns
df <- data.frame(matrix(ncol = 6, nrow = 7*30))
colnames(df) <- columns
organ_size <- sum(abundance_data_for_estimation$Count)
for(k in 1:10) {
  for (j in 1:length(sample_fraction)) {
    freq_df <- data.frame(Value = integer(0)) #This collects ABUNDANCE data
    for (i in 1:n_samples) {
      sampled_cells <- sample(abundance_data_for_estimation$`Cell type`, size = organ_size*sample_fraction[j], replace = TRUE, prob = abundance_data_for_estimation$prob)
      sample_freq_df <- generate_freq_table(sampled_cells, i)
      freq_df <- merge(freq_df, sample_freq_df, by = "Value", all = TRUE)
      freq_df[is.na(freq_df)] <- 0
    }
    #freq_df <- freq_df[freq_df$Frequency_Sample_1 != 0, ]
    freq_df <- freq_df[,-1]
    freq_df <- freq_df[rowSums(freq_df == 0) != ncol(freq_df), ]
    
    freq_df <- freq_df %>%
      mutate(across(everything(), ~ ifelse(. > 0, 1, .)))
    
    row_index <- (k - 1) * length(sample_fraction) + j
    # Add estimates from estimates_df to new_row
    df[row_index, 1] <- sample_fraction[j]*100
    df[row_index,"Found in sample"] <- nrow(freq_df)
    df[row_index,"ICE"] <- ICE(freq_df)
    df[row_index,"Chao2"] <- chao2(freq_df)
    df[row_index,"Jack1"] <- jack1(freq_df, abund = FALSE)
    df[row_index,"Jack2"] <- jack2(freq_df, abund = FALSE)
  }
}
df[ ,"Ground truth"] <- nrow(abundance_data_for_estimation)

df_long <- pivot_longer(df, 
                        cols = -c(`Sample coverage`, `Ground truth`), 
                        names_to = "Estimate_Type", 
                        values_to = "Estimate_Value")

ggplot(df_long, aes(x = `Sample coverage`, y = Estimate_Value, color = Estimate_Type)) +
  geom_smooth(method = "loess", se = TRUE, size = 2) +  # Smooth the lines with LOESS
  geom_hline(aes(yintercept = `Ground truth`), color = "red", linetype = "dashed") +
  scale_color_brewer(palette = "Paired") +
  scale_x_log10() +  # Log scale for x-axis if desired
  labs(title = "Estimates vs Sample Fraction",
       subtitle = "Comparison of Different Estimation Methods",
       x = "Sample coverage %",
       y = "Estimated cellular richness",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.major = element_line(color = "gray90"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 12, colour = "black"),  # Consistent y-axis label size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"), # Rotate and style x-axis labels
        axis.title.x = element_text(size = 16),                  # Standardize x-axis title size
        axis.title.y = element_text(size = 16),                  # Standardize y-axis title size
        plot.title = element_text(size = 14, face = "bold"))

ggsave("Largest_organ_estimations_all.svg", scale = 1.5, device = "svg")

df_long[ ,"Error"] <- (df_long[ ,"Estimate_Value"] - df_long[ ,"Ground truth"]) / df_long[ ,"Ground truth"] * 100

ggplot(df_long, aes(x = `Sample coverage`, y = Error, color = Estimate_Type)) +
  geom_smooth(method = "loess", se = TRUE, size = 2) +  # Smooth the lines with LOESS
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  scale_color_brewer(palette = "Paired") +
  scale_x_log10() +  # Log scale for x-axis if desired
  labs(title = "Estimation error",
       subtitle = "Comparison of Different Estimation Methods",
       x = "Sample coverage %",
       y = "Error %",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.major = element_line(color = "gray90"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 12, colour = "black"),  # Consistent y-axis label size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"), # Rotate and style x-axis labels
        axis.title.x = element_text(size = 16),                  # Standardize x-axis title size
        axis.title.y = element_text(size = 16),                  # Standardize y-axis title size
        plot.title = element_text(size = 14, face = "bold"))

ggsave("Largest_organ_estimations_error_all.svg", scale = 1.5, device = "svg")

ggplot(df_long, aes(x = Estimate_Type, y = Error, color = Estimate_Type)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 4) +  # Add individual points (dots) with jitter
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +  # Mark median
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +  # Reference line at 0
  scale_color_brewer(palette = "Paired") +  # Use the "Paired" color palette
  labs(title = "Comparison of Error Measurements Across Estimation Methods",
       x = "Estimation Method",
       y = "Error (%)",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis labels
        plot.title = element_text(size = 16, face = "bold"))

ggsave("Largest_organ_all_error.svg", scale = 1, device = "svg")


#Estimation for 3-organs
oneorgan <- as.data.frame(scaled_dataframe$`20`)
colnames(oneorgan) <- "Count"
oneorgan <- oneorgan[oneorgan$Count != 0, ]
oneorgan <- as.data.frame(oneorgan)
colnames(oneorgan) <- "Count"
oneorgan <- round(oneorgan)
oneorgan$"Cell type" <- 1:nrow(oneorgan)
oneorgan$"prob" <- oneorgan$Count / sum(oneorgan$Count)
abundance_data_for_estimation <- oneorgan

#Abundance dataframe of the population
abundance_data_for_estimation$Count <- round(abundance_data_for_estimation$Count)
abundance_data_for_estimation$prob <- abundance_data_for_estimation$Count / sum(abundance_data_for_estimation$Count)

generate_freq_table <- function(sample_data, sample_num) {
  observed_freq <- table(sample_data)
  observed_freq_df <- as.data.frame(observed_freq)
  
  # Ensure columns are named correctly
  colnames(observed_freq_df) <- c("Value", paste0("Frequency_Sample_", sample_num))
  
  return(observed_freq_df)
}

# Define the column names
columns <- c("Sample coverage", "Found in sample", "ICE", "Chao2", "Jack1", "Jack2")

# Create an empty data frame with these columns
df <- data.frame(matrix(ncol = 6, nrow = 7*30))
colnames(df) <- columns
for(k in 1:10) {
  for (j in 1:length(sample_fraction)) {
    freq_df <- data.frame(Value = integer(0)) #This collects ABUNDANCE data
    for (i in 1:n_samples) {
      sampled_cells <- sample(abundance_data_for_estimation$`Cell type`, size = organ_size*sample_fraction[j], replace = TRUE, prob = abundance_data_for_estimation$prob)
      sample_freq_df <- generate_freq_table(sampled_cells, i)
      freq_df <- merge(freq_df, sample_freq_df, by = "Value", all = TRUE)
      freq_df[is.na(freq_df)] <- 0
    }
    #freq_df <- freq_df[freq_df$Frequency_Sample_1 != 0, ]
    freq_df <- freq_df[,-1]
    freq_df <- freq_df[rowSums(freq_df == 0) != ncol(freq_df), ]
    
    freq_df <- freq_df %>%
      mutate(across(everything(), ~ ifelse(. > 0, 1, .)))
    
    row_index <- (k - 1) * length(sample_fraction) + j
    # Add estimates from estimates_df to new_row
    df[row_index, 1] <- sample_fraction[j]*100
    df[row_index,"Found in sample"] <- nrow(freq_df)
    df[row_index,"ICE"] <- ICE(freq_df)
    df[row_index,"Chao2"] <- chao2(freq_df)
    df[row_index,"Jack1"] <- jack1(freq_df, abund = FALSE)
    df[row_index,"Jack2"] <- jack2(freq_df, abund = FALSE)
  }
}
df[ ,"Ground truth"] <- nrow(abundance_data_for_estimation)

df_long <- pivot_longer(df, 
                        cols = -c(`Sample coverage`, `Ground truth`), 
                        names_to = "Estimate_Type", 
                        values_to = "Estimate_Value")

ggplot(df_long, aes(x = `Sample coverage`, y = Estimate_Value, color = Estimate_Type)) +
  geom_smooth(method = "loess", se = TRUE, size = 2) +  # Smooth the lines with LOESS
  geom_hline(aes(yintercept = `Ground truth`), color = "red", linetype = "dashed") +
  scale_color_brewer(palette = "Paired") +
  scale_x_log10() +  # Log scale for x-axis if desired
  labs(title = "Estimates vs Sample Fraction",
       subtitle = "Comparison of Different Estimation Methods",
       x = "Sample coverage %",
       y = "Estimated cellular richness",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.major = element_line(color = "gray90"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 12, colour = "black"),  # Consistent y-axis label size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"), # Rotate and style x-axis labels
        axis.title.x = element_text(size = 16),                  # Standardize x-axis title size
        axis.title.y = element_text(size = 16),                  # Standardize y-axis title size
        plot.title = element_text(size = 14, face = "bold"))

ggsave("Medium_organ_estimations_all.svg", scale = 1.5, device = "svg")

df_long[ ,"Error"] <- (df_long[ ,"Estimate_Value"] - df_long[ ,"Ground truth"]) / df_long[ ,"Ground truth"] * 100

ggplot(df_long, aes(x = `Sample coverage`, y = Error, color = Estimate_Type)) +
  geom_smooth(method = "loess", se = TRUE, size = 2) +  # Smooth the lines with LOESS
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  scale_color_brewer(palette = "Paired") +
  scale_x_log10() +  # Log scale for x-axis if desired
  labs(title = "Estimation error",
       subtitle = "Comparison of Different Estimation Methods",
       x = "Sample coverage %",
       y = "Error %",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.major = element_line(color = "gray90"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 12, colour = "black"),  # Consistent y-axis label size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"), # Rotate and style x-axis labels
        axis.title.x = element_text(size = 16),                  # Standardize x-axis title size
        axis.title.y = element_text(size = 16),                  # Standardize y-axis title size
        plot.title = element_text(size = 14, face = "bold"))

ggsave("Medium_organ_estimations_error_all.svg", scale = 1.5, device = "svg")

ggplot(df_long, aes(x = Estimate_Type, y = Error, color = Estimate_Type)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +  # Add individual points (dots) with jitter
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +  # Mark median
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +  # Reference line at 0
  scale_color_brewer(palette = "Paired") +  # Use the "Paired" color palette
  labs(title = "Comparison of Error Measurements Across Estimation Methods",
       x = "Estimation Method",
       y = "Error (%)",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis labels
        plot.title = element_text(size = 16, face = "bold"))

ggsave("Smallest_organ_all_error.svg", scale = 1.5, device = "svg")

#Estimation for 3-organs
oneorgan <- as.data.frame(scaled_dataframe$`2`)
colnames(oneorgan) <- "Count"
oneorgan <- oneorgan[oneorgan$Count != 0, ]
oneorgan <- as.data.frame(oneorgan)
colnames(oneorgan) <- "Count"
oneorgan <- round(oneorgan)
oneorgan$"Cell type" <- 1:nrow(oneorgan)
oneorgan$"prob" <- oneorgan$Count / sum(oneorgan$Count)
abundance_data_for_estimation <- oneorgan

#Abundance dataframe of the population
abundance_data_for_estimation$Count <- round(abundance_data_for_estimation$Count)
abundance_data_for_estimation$prob <- abundance_data_for_estimation$Count / sum(abundance_data_for_estimation$Count)

generate_freq_table <- function(sample_data, sample_num) {
  observed_freq <- table(sample_data)
  observed_freq_df <- as.data.frame(observed_freq)
  
  # Ensure columns are named correctly
  colnames(observed_freq_df) <- c("Value", paste0("Frequency_Sample_", sample_num))
  
  return(observed_freq_df)
}

# Define the column names
columns <- c("Sample coverage", "Found in sample", "ICE", "Chao2", "Jack1", "Jack2")

# Create an empty data frame with these columns
df <- data.frame(matrix(ncol = 6, nrow = 7*30))
colnames(df) <- columns
for(k in 1:10) {
  for (j in 1:length(sample_fraction)) {
    freq_df <- data.frame(Value = integer(0)) #This collects ABUNDANCE data
    for (i in 1:n_samples) {
      sampled_cells <- sample(abundance_data_for_estimation$`Cell type`, size = organ_size*sample_fraction[j], replace = TRUE, prob = abundance_data_for_estimation$prob)
      sample_freq_df <- generate_freq_table(sampled_cells, i)
      freq_df <- merge(freq_df, sample_freq_df, by = "Value", all = TRUE)
      freq_df[is.na(freq_df)] <- 0
    }
    #freq_df <- freq_df[freq_df$Frequency_Sample_1 != 0, ]
    freq_df <- freq_df[,-1]
    freq_df <- freq_df[rowSums(freq_df == 0) != ncol(freq_df), ]
    
    freq_df <- freq_df %>%
      mutate(across(everything(), ~ ifelse(. > 0, 1, .)))
    
    row_index <- (k - 1) * length(sample_fraction) + j
    # Add estimates from estimates_df to new_row
    df[row_index, 1] <- sample_fraction[j]*100
    df[row_index,"Found in sample"] <- nrow(freq_df)
    df[row_index,"ICE"] <- ICE(freq_df)
    df[row_index,"Chao2"] <- chao2(freq_df)
    df[row_index,"Jack1"] <- jack1(freq_df, abund = FALSE)
    df[row_index,"Jack2"] <- jack2(freq_df, abund = FALSE)
  }
}
df[ ,"Ground truth"] <- nrow(abundance_data_for_estimation)

df_long <- pivot_longer(df, 
                        cols = -c(`Sample coverage`, `Ground truth`), 
                        names_to = "Estimate_Type", 
                        values_to = "Estimate_Value")

ggplot(df_long, aes(x = `Sample coverage`, y = Estimate_Value, color = Estimate_Type)) +
  geom_smooth(method = "loess", se = TRUE, size = 2) +  # Smooth the lines with LOESS
  geom_hline(aes(yintercept = `Ground truth`), color = "red", linetype = "dashed") +
  scale_color_brewer(palette = "Paired") +
  scale_x_log10() +  # Log scale for x-axis if desired
  labs(title = "Estimates vs Sample Fraction",
       subtitle = "Comparison of Different Estimation Methods",
       x = "Sample coverage %",
       y = "Estimated cellular richness",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.major = element_line(color = "gray90"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 12, colour = "black"),  # Consistent y-axis label size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"), # Rotate and style x-axis labels
        axis.title.x = element_text(size = 16),                  # Standardize x-axis title size
        axis.title.y = element_text(size = 16),                  # Standardize y-axis title size
        plot.title = element_text(size = 14, face = "bold"))

ggsave("Medium_organ_estimations_all.svg", scale = 1.5, device = "svg")

df_long[ ,"Error"] <- (df_long[ ,"Estimate_Value"] - df_long[ ,"Ground truth"]) / df_long[ ,"Ground truth"] * 100

ggplot(df_long, aes(x = `Sample coverage`, y = Error, color = Estimate_Type)) +
  geom_smooth(method = "loess", se = TRUE, size = 2) +  # Smooth the lines with LOESS
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  scale_color_brewer(palette = "Paired") +
  scale_x_log10() +  # Log scale for x-axis if desired
  labs(title = "Estimation error",
       subtitle = "Comparison of Different Estimation Methods",
       x = "Sample coverage %",
       y = "Error %",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.major = element_line(color = "gray90"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 12, colour = "black"),  # Consistent y-axis label size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"), # Rotate and style x-axis labels
        axis.title.x = element_text(size = 16),                  # Standardize x-axis title size
        axis.title.y = element_text(size = 16),                  # Standardize y-axis title size
        plot.title = element_text(size = 14, face = "bold"))     # Bold title

ggsave("Medium_organ_estimations_error_all.svg", scale = 1.5, device = "svg")

ggplot(df_long, aes(x = Estimate_Type, y = Error, color = Estimate_Type)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +  # Add individual points (dots) with jitter
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +  # Mark median
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +  # Reference line at 0
  scale_color_brewer(palette = "Paired") +  # Use the "Paired" color palette
  labs(title = "Comparison of Error Measurements Across Estimation Methods",
       x = "Estimation Method",
       y = "Error (%)",
       color = "Estimation Method") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis labels
        plot.title = element_text(size = 16, face = "bold"))

ggsave("Medium_organ_all_error.svg", scale = 1.5, device = "svg")