############################################################
## Faceted rank–abundance barplots (NO p-value, NO fit line)
## X = cell types ranked by decreasing abundance (rank)
## Y = abundance (count)
############################################################

pkgs <- c("openxlsx","dplyr","purrr","tibble","ggplot2","stringr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install)) install.packages(to_install)

library(openxlsx)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(stringr)

## ---- File paths (EDIT if needed) ----
file_paths <- list(
  Eye = "TS eye.xlsx",
  Bladder = "TS bladder.xlsx",
  Heart = "TS heart.xlsx",
  Kidney = "TS kidney.xlsx",
  "Bone marrow" = "TS bone marrow.xlsx",
  Liver = "TS liver.xlsx",
  "Salivary gland" = "TS salivary gland.xlsx",
  Lung = "TS lung.xlsx",
  "Lymph node" = "TS lymph node.xlsx",
  "Mammary gland" = "TS mammary gland.xlsx",
  Pancreas = "TS pancreas.xlsx",
  Prostate = "TS prostate.xlsx",
  Skin = "TS skin.xlsx",
  LargeIntestine = "TS large intestine.xlsx",
  SmallIntestine = "TS small intestine.xlsx",
  Spleen = "TS spleen.xlsx",
  Thymus = "TS thymus.xlsx",
  Tongue = "TS tongue.xlsx",
  Trachea = "TS trachea.xlsx",
  Uterus = "TS uterus.xlsx"
)

missing <- file_paths[!file.exists(unlist(file_paths))]
if(length(missing) > 0){
  stop("These files do not exist:\n", paste(unname(missing), collapse="\n"))
}

## ---- Helpers: robust column detection ----
norm_name <- function(x){
  x <- as.character(x)
  x2 <- iconv(x, to="ASCII//TRANSLIT")
  x2[is.na(x2)] <- x[is.na(x2)]
  x2 <- tolower(x2)
  gsub("[^a-z0-9]", "", x2)
}

get_required_cols <- function(df, file_name){
  nms <- names(df)
  nn  <- norm_name(nms)

  ct_hit <- which(str_detect(nn, "celltype") | str_detect(nn, "elltype"))
  c_hit  <- which(str_detect(nn, "count") | str_detect(nn, "ount") | nn == "n")

  if(length(ct_hit) < 1 || length(c_hit) < 1){
    stop(
      "File ", file_name, " must contain a cell type column and a count column.\n",
      "Found: ", paste(nms, collapse=", ")
    )
  }
  list(ct = nms[ct_hit[1]], count = nms[c_hit[1]])
}

## ---- Load all data ----
dat <- imap_dfr(file_paths, function(fp, organ_name){
  df <- read.xlsx(fp)
  cols <- get_required_cols(df, basename(fp))

  raw_count <- df[[cols$count]]
  raw_count <- if(is.factor(raw_count)) as.character(raw_count) else raw_count
  raw_count <- gsub(",", "", as.character(raw_count))
  count_num <- suppressWarnings(as.numeric(raw_count))

  tibble(
    organ = organ_name,
    cell_type = df[[cols$ct]],
    count = count_num
  ) %>%
    filter(is.finite(count), count >= 0) %>%
    mutate(count = as.integer(round(count)))
})

## ---- Create rank per organ ----
rank_df <- dat %>%
  group_by(organ) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

## ---- Plot: bars only ----
p <- ggplot(rank_df, aes(x = rank, y = count)) +
  geom_col(width = 0.9, fill = "grey40") +
  facet_wrap(~ organ, scales = "free") +
  labs(
    title = "Cell-type abundance distributions",
    x = "Cell types (ranked by decreasing abundance)",
    y = "Abundance"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey80", colour = "grey40"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("celltype_rank_abundance_bars_only.png", p, width = 14, height = 8, dpi = 300)
ggsave("celltype_rank_abundance_bars_only.pdf", p, width = 14, height = 8)

p
