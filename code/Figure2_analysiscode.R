#Combine TS data from 20 organs
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
options(scipen=10)


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

#How are non-diagonal values distributed?
upper_tri_values <- correlation_matrix[upper.tri(correlation_matrix, diag = FALSE)]

# Print the resulting vector
print(upper_tri_values)

upper_tri_df <- data.frame(Value = upper_tri_values)

loglinear_model <- glm(Value ~ 1, family = gaussian(link = "log"), data = upper_tri_df)
summary(loglinear_model)
# Extracting the meanlog (the coefficient of the intercept)
meanlog <- coef(loglinear_model)[1]

# Calculating the standard deviation of the log-transformed residuals (sdlog)
# Residual standard error can be extracted from the model summary
model_summary <- summary(loglinear_model)
sdlog <- model_summary$dispersion^0.5  # Residual standard deviation

# Create a density plot
ggplot(as.data.frame(upper_tri_df), aes(x = Value)) +
  geom_density(fill = "blue", alpha = 0.7) +
  labs(title = "Density Plot of Upper Triangular Values in Correlation Matrix",
       x = "Correlation Value",
       y = "Density") +
  theme_minimal()

ggsave("Nondiagonal_correlation_distribution.svg", scale = 1.5, device = "svg")

#log-transformed data
epsilon <- 1e-5
data_matrix[data_matrix == 0] <- epsilon
log_data_matrix <- log(data_matrix)

#Cov matrix
log_cov_matrix <- cov(log_data_matrix)

#How log-cov_matrix is distributed
upper_tri_values <- log_cov_matrix[upper.tri(log_cov_matrix, diag = FALSE)]
diagonal_values <- diag(log_cov_matrix)
upper_tri_df <- data.frame(Value = upper_tri_values)

# Create a density plot
ggplot(as.data.frame(upper_tri_df), aes(x = Value)) +
  geom_density(fill = "blue", alpha = 0.7) +
  labs(title = "Density Plot of Upper Triangular Values in COV Matrix",
       x = "covariance Values",
       y = "Density") +
  theme_minimal()

lm_model <- lm(upper_tri_values ~ 1)

# Extracting the mean and standard deviation from the model
mean_value <- coef(lm_model)  # This gives the mean
residuals_sd <- summary(lm_model)$sigma  # This gives the standard deviation (of residuals)

mean_value
residuals_sd

# Summary of the fitted model
summary(lm_model)

# Compute the mean vector and covariance matrix
log_means <- colMeans(log_data_matrix)
sim_log_means <- rnorm(20, -9, 0.84)

#Building sim_correlation_matrix (other page)
num_organs <- 20  # Dimension of the covariance matrix
mean <- 6.448077   # Mean parameter for the lognormal distribution
log <- 4.98086  # Standard deviation parameter for the lognormal distribution
log_cov_matrix_20 <- generate_cov_matrix(num_organs, mean, log)
print(log_cov_matrix_20)

# Simulate log-normal data for 500 cell types
simulated_normal <- mvrnorm(n=700, mu=sim_log_means, Sigma=log_cov_matrix_20)
simulated_lognormal <- exp(simulated_normal)


# Scale the simulated data to match the original data's row sums
#original_row_sums <- rowSums(data_matrix)
#scaled_simulated_data <- t(apply(simulated_lognormal, 1, function(row) row / sum(row) * mean(original_row_sums)))

# Convert the scaled data to integers and set very small values to zero
integer_simulated_data <- round(simulated_lognormal)
integer_simulated_data[integer_simulated_data < epsilon] <- 0

total_cells_in_20_organs <- 10000000000000

scaling_factor <- total_cells_in_20_organs / sum(colSums(integer_simulated_data))
scaled_dataframe <- integer_simulated_data * scaling_factor

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

#Save figure 2C
ggsave("sim_human_unique_shared.svg", scale = 1.5, device = "svg")

# Calculate the correlation matrix
correlation_matrix <- cor(as.matrix(scaled_dataframe), use = "pairwise.complete.obs")


#Save figure 2D
corrplot::corrplot(correlation_matrix, method = "color", type = "upper", xlab = "testx", ylab = "testy", tl.col = "black", tl.srt = 45)
ggsave("Organ_correlation.svg", scale = 1.5, device = "svg")

#This data has cell types
nrow(abundance_data)

#Estimation for 3-organs
oneorgan <- as.data.frame(scaled_dataframe$`10`)
colnames(oneorgan) <- "Count"
oneorgan <- oneorgan[oneorgan$Count != 0, ]
oneorgan <- as.data.frame(oneorgan)
colnames(oneorgan) <- "Count"
oneorgan <- round(oneorgan)
oneorgan$"Cell type" <- 1:nrow(oneorgan)
oneorgan$"prob" <- oneorgan$Count / sum(oneorgan$Count)
abundance_data_for_estimation <- oneorgan

n_samples <- 6
sample_fraction <- c(1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05)

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

ggsave("Largest_organ_all_error.svg", scale = 1.5, device = "svg")


#Estimation for 3-organs
oneorgan <- as.data.frame(scaled_dataframe$`3`)
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

ggsave("Smallest_organ_estimations_all.svg", scale = 1.5, device = "svg")

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

ggsave("Smallest_organ_estimations_error_all.svg", scale = 1.5, device = "svg")

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

