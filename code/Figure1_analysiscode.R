#Libraries
library(UpSetR)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(knitr)
library(kableExtra)
library(gridExtra)
library(webshot2)
library(stats)
library(fitdistrplus)
library(openxlsx)
library(tidyr)
library(readxl)
library(reshape2)


# Define the file paths for the Excel files
file_paths <- list(
  Eye = "Eye - Tabula Sapiens.xlsx",
  Bladder = "Bladder - Tabula Sapiens.xlsx",
  Heart = "Heart - Tabula Sapiens.xlsx",
  Kidney = "Kidney - Tabula Sapiens.xlsx",
  "Large intestine" = "Large intestine - Tabula Sapiens.xlsx",
  "Bone marrow" = "Bone marrow - Tabula Sapiens.xlsx",
  Liver = "Liver - Tabula Sapiens.xlsx",
  "Salivary gland" = "Salivary gland - Tabula Sapiens.xlsx",
  Lung = "Lung - Tabula Sapiens.xlsx",
  "Lymph node" = "Lymph node - Tabula Sapiens.xlsx",
  "Mammary gland" = "Mammary gland - Tabula Sapiens.xlsx",
  Pancreas = "Pancreas - Tabula Sapiens.xlsx",
  Prostate = "Prostate - Tabula Sapiens.xlsx",
  Skin = "Skin - Tabula Sapiens.xlsx",
  "Small intestine" = "Small intestine - Tabula Sapiens.xlsx",
  Spleen = "Spleen - Tabula Sapiens.xlsx",
  Thymus = "Thymus - Tabula Sapiens.xlsx",
  Tongue = "Tongue - Tabula Sapiens.xlsx",
  Trachea = "Trachea - Tabula Sapiens.xlsx",
  Uterus = "Uterus - Tabula Sapiens.xlsx"
  
  
)

# Check if files exist
for (file_path in file_paths) {
  if (!file.exists(file_path)) {
    stop(paste("File does not exist:", file_path))
  }
}

# Read the data from each Excel file and extract the cell types
cell_data <- lapply(file_paths, function(file_path) {
  data <- read.xlsx(file_path)
  cell_types <- data[[1]]  #first column contains the cell types
  return(cell_types)
})

# Convert the list of cell types into a binary matrix
all_cell_types <- unique(unlist(cell_data))
binary_matrix <- sapply(cell_data, function(cell_types) {
  as.integer(all_cell_types %in% cell_types)
})
colnames(binary_matrix) <- names(file_paths)
binary_matrix <- as.data.frame(binary_matrix)
binary_matrix$cell_type <- all_cell_types


# Create the UpSet plot
# Set the file name and path for the SVG file
svg_filename <- "UpSet_plot.svg"

# Open the SVG device
svg(svg_filename, width = 10, height = 7)

# Create the UpSet plot
upset(binary_matrix, 
      sets = names(file_paths), 
      order.by = "freq")

# Close the SVG device
dev.off()

#Figure 1A
upset(binary_matrix, sets = names(file_paths), order.by = "freq")

# Ensure binary_matrix is a data frame
binary_matrix <- as.data.frame(binary_matrix)

# Calculate total cell types per organ
organ_cell_counts <- binary_matrix %>%
  dplyr::select(-cell_type) %>%  # Remove the cell_type column
  summarise(across(everything(), ~ sum(. == 1)))  # Calculate total cell types for each organ

# Transpose the data and convert to data frame
organ_cell_counts <- as.data.frame(t(organ_cell_counts))
colnames(organ_cell_counts) <- "total"
organ_cell_counts$organ <- rownames(organ_cell_counts)

# Calculate unique cell types per organ
unique_counts <- sapply(binary_matrix %>% dplyr::select(-cell_type), function(col) {
  sum(col == 1 & rowSums(binary_matrix %>% dplyr::select(-cell_type)) == 1)
})
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




ggsave("Cell type uniqueness stacked barchart.svg", scale = 1.5, device = "svg")

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


corrplot::corrplot(correlation_matrix, method = "color", type = "upper", xlab = "testx", ylab = "testy", tl.col = "black", tl.srt = 45)
ggsave("Organ_correlation.svg", scale = 1.5, device = "svg")
  
