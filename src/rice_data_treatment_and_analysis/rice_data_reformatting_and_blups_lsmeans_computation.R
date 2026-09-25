# script meant to reformat data for genomic prediction
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(devtools)
install_other_requirements <- F
if (install_other_requirements) {
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("snpStats")
  BiocManager::install("mixOmicsTeam/mixOmics")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
}
library(mixOmics)
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(Matrix)
library(graphics)
library(htmlwidgets)
library(rstudioapi)
library(stringr)
library(tidyr)
library(dplyr)
library(lsmeans)
library(lme4)
library(pheatmap)

# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../functions.R")

# set options to increase memory
options(expressions = 5e5)
options(warn = -1)
emm_options(rg.limit = 10e6)

# set path for genomic data and phenotype data
genom_dir_path <- "../../data/genomic_data/"
pheno_dir_path <- "../../data/phenotype_data/"

# set output result path for genomic graphics
output_genom_graphics_path <- "../../results/genomic_prediction_graphics/"

# set maximum number of principal components to be tested using akaike
# information criterion
max_n_comp_ <- 10

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 15
min_dist_ <- 0.1

# define traits_
traits_ <- c("FL", "PH", "YLD", "ZN")

# get genomic data
geno_df <- t(as.data.frame(fread(paste0(
  genom_dir_path,
  "PCT27_TP334_genotypes.txt"
))))
col_names_ <- geno_df[1, -1]
row_names_ <- geno_df[-1, 1]
geno_df <- geno_df[-1, -1]
colnames(geno_df) <- col_names_
rownames(geno_df) <- row_names_

# get phenotype data
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "PCT27_TP334_phenotypes.csv"
)))
colnames(pheno_df)[match(
  "DNAID",
  colnames(pheno_df)
)] <- "Genotype"

# get common genotypes between genomic and phenotype data
pheno_df <- match_indices(pheno_df, geno_df)
geno_df <- geno_df[rownames(geno_df) %in% pheno_df$Genotype, ]

# create new factor generation-block-trial
pheno_df$Envir <- paste0(
  pheno_df$TRIAL, "_",
  pheno_df$BLOC, "_",
  pheno_df$GENERATION
)

# write reformatted datasets
fwrite(pheno_df,
  file = paste0(pheno_dir_path, "phenotype_data.csv")
)
fwrite(as.data.frame(geno_df),
  file = paste0(genom_dir_path, "genomic_data.csv"),
  row.names = T
)

# convert geno_df to numeric matrix
geno_mat <- scale(apply(geno_df, 2, as.numeric), center = T, scale = F)
k_mat <- crossprod(t(geno_mat))

# create a heatmap for genomic covariance matrix
png(
  paste0(
    output_genom_graphics_path,
    "rice_genomic_covariance_heatmap.png"
  ),
  width = 1200, height = 1200, res = 150
)
pheatmap(k_mat,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Rice dataset genomic covariance matrix heatmap",
  color = colorRampPalette(c("blue", "red"))(100),
  show_rownames = T,
  show_colnames = T
)
dev.off()

# compute pca for genomic data
geno_pca <- mixOmics::pca(apply(geno_df, 2, as.numeric), 2,
  center = T,
  scale = F,
  ncomp = max_n_comp_
)
pc_coord_df_ <- as.data.frame(geno_pca$variates)[, 1:max_n_comp_]
colnames(pc_coord_df_) <- str_replace_all(
  colnames(pc_coord_df_),
  pattern = "X.", replacement = ""
)
pc_var_names_ <- colnames(pc_coord_df_)
pc_coord_df_$Genotype <- rownames(geno_df)
geno_pca_exp_var_ <- geno_pca$explained_variance

# create pca plot
pca_fig_x_y <- plot_ly(
  data = pc_coord_df_,
  x = ~PC1, y = ~PC2,
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for rice dataset genomic data (9928 SNP)",
    xaxis = list(title = paste0(
      names(geno_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(geno_pca_exp_var_)[1], 2), "%"
    )),
    yaxis = list(title = paste0(
      names(geno_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(geno_pca_exp_var_)[2], 2), "%"
    ))
  )

# save pca graphic
saveWidget(pca_fig_x_y, file = paste0(
  output_genom_graphics_path, "rice_pca_genomic_data.html"
))

# compute umap in 2d
geno_umap_2d_model <- umap(
  apply(geno_df, 2, as.numeric),
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)
geno_umap_2d <- data.frame(geno_umap_2d_model[["layout"]])

# create umap plot
umap_fig_x_y <- plot_ly(
  data = geno_umap_2d,
  x = ~X1, y = ~X2,
  type = "scatter", mode = "markers", color = "orange"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "UMAP 2D plot for rice dataset genomic data (9928 SNP)",
    xaxis = list(title = "First component"),
    yaxis = list(title = "Second component")
  )

# save umap graphic
saveWidget(umap_fig_x_y, file = paste0(
  output_genom_graphics_path, "rice_umap_genomic_data.html"
))

# get ls-means and blups (using pca as fixed covariates for population structure
# correction) for each trait across all generations, blocks and trials

# initialize list for traits ls-means and blups
blup_list_ <- vector("list", length(traits_))
blup_pca_list_ <- vector("list", length(traits_))
ls_means_list_ <- vector("list", length(traits_))
names(blup_list_) <- traits_
names(blup_pca_list_) <- traits_
names(ls_means_list_) <- traits_
aic_ <- rep(0, max_n_comp_)

# compute ls-means and blups
for (trait_ in traits_) {
  pheno_df_trait_ <- pheno_df[
    which(pheno_df[, trait_] != "na"),
  ]
  pheno_df_trait_[, trait_] <- as.numeric(pheno_df_trait_[, trait_])

  # merge trait individual phenotypes with associated
  # pc coordinates for genotypes
  pheno_df_trait_ <- merge(pheno_df_trait_, pc_coord_df_,
    by = "Genotype", all = TRUE
  )

  # compute aic values in order to select number of pcs
  for (n_comp_ in 1:max_n_comp_) {
    lmer_pca_model_ <- lmer(
      as.formula(paste0(
        trait_,
        " ~ 1 + Envir + ", paste(pc_var_names_[1:n_comp_],
          collapse = " + "
        ),
        " + (1 | Genotype)"
      )),
      data = pheno_df_trait_,
      REML = FALSE
    )
    aic_[n_comp_] <- AIC(lmer_pca_model_)
  }
  n_opt_comp_aic_ <- which.min(aic_)
  print(paste0("number of pc selected: ", n_opt_comp_aic_))

  # estimate model based on selected number of pcs which minimize aic
  lmer_pca_model_ <- lmer(
    as.formula(paste0(
      trait_,
      " ~ 1 + Envir + ", paste(pc_var_names_[1:n_opt_comp_aic_],
        collapse = " + "
      ),
      " + (1 | Genotype)"
    )),
    data = pheno_df_trait_,
    REML = TRUE
  )
  df_ <- data.frame(
    "Genotype" = rownames(ranef(lmer_pca_model_)$Genotype),
    "blup_pca" = as.numeric(unlist(ranef(lmer_pca_model_)$Genotype))
  )
  blup_pca_list_[[trait_]] <- df_

  # estimate model based on linear mixed model (LMM)
  lmer_model_ <- lmer(
    as.formula(paste0(
      trait_,
      " ~ 1 + Envir + (1 | Genotype)"
    )),
    data = pheno_df_trait_,
    REML = TRUE
  )
  df_ <- data.frame(
    "Genotype" = rownames(ranef(lmer_model_)$Genotype),
    "blup" = as.numeric(unlist(ranef(lmer_model_)$Genotype))
  )
  blup_list_[[trait_]] <- df_

  # apply multiple regression first
  lm_ <- lm(
    as.formula(
      paste0(
        trait_,
        " ~  1 + Genotype + Envir"
      )
    ),
    data = pheno_df_trait_
  )
  # get genotype ls-means for trait_
  ls_means_list_[[trait_]] <- as.data.frame(lsmeans(
    lm_,
    ~Genotype
  ))[, c("Genotype", "lsmean")]
}

# reduce lsmeans list
ls_means_df <- Reduce(
  function(x, y) {
    merge(x, y,
      by = "Genotype",
      all = T
    )
  },
  ls_means_list_
)
colnames(ls_means_df) <- c("Genotype", traits_)

# reduce blup list
blup_pca_df <- Reduce(
  function(x, y) {
    merge(x, y,
      by = "Genotype",
      all = T
    )
  },
  blup_pca_list_
)
colnames(blup_pca_df) <- c("Genotype", traits_)

# reduce blup list
blup_df <- Reduce(
  function(x, y) {
    merge(x, y,
      by = "Genotype",
      all = T
    )
  },
  blup_list_
)
colnames(blup_df) <- c("Genotype", traits_)

# write ls-means
fwrite(ls_means_df, file = paste0(
  pheno_dir_path,
  "ls_mean_breeding_values.csv"
))

# write blups
fwrite(blup_pca_df, file = paste0(
  pheno_dir_path,
  "blup_pca_breeding_values.csv"
))

# write blups
fwrite(blup_df, file = paste0(
  pheno_dir_path,
  "blup_breeding_values.csv"
))
