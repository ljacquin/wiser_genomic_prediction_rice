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
pheno_df$GENERATION_BLOC_TRIAL <- paste0(
  pheno_df$GENERATION, "_",
  pheno_df$BLOC,
  "_", pheno_df$TRIAL
)

# write reformatted datasets
fwrite(pheno_df,
       file = paste0(pheno_dir_path, "phenotype_data.csv")
)
fwrite(as.data.frame(geno_df),
       file = paste0(genom_dir_path, "genomic_data.csv"),
       row.names = T
)
