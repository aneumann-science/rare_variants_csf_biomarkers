# Load libraries
library(SKAT)
library(MetaSKAT)
library(data.table)
library(dplyr)
library(pbmcapply)

# Load ADNI phenotype and covariate data
load("phenotypes/adni/mmse_adni.Rdata")
load("phenotypes/adni/pca_adni.Rdata")
load("phenotypes/adni/csf_cov_diag_adni.Rdata")

# Load EMIF phenotype and covariate data
load("phenotypes/emif/mmse_emif.Rdata")
load("phenotypes/emif/pca_emif.Rdata")
load("phenotypes/emif/csf_cov_diag_emif.Rdata")

# Read in paths for available genes
genes_adni <- list.files(path="../adni_emif_genotypes/adni/protein_coding/plink", pattern="*.raw", full.names=FALSE, recursive=FALSE)
genes_emif <- list.files(path="../adni_emif_genotypes/emif/protein_coding/plink", pattern="*.raw", full.names=FALSE, recursive=FALSE)

# Match genes present in both
genes_adni <- genes_adni[genes_adni %in% genes_emif]
genes_emif <- genes_emif[genes_emif %in% genes_adni]

# Load SKAT-O Meta function
source("skat_o_meta_diag.R")

####### PCA
# Tau
tau_pc_meta_p_values.list <- pbmclapply(genes_adni, function(gene) {
  skat_o_meta(pca_adni.data,"RC1",pca_emif.data,"RC1",csf_cov_diag_adni.data,csf_cov_diag_emif.data,gene)
}, mc.cores = 12)

tau_pc_meta_results.data <- do.call(rbind, tau_pc_meta_p_values.list)
save(tau_pc_meta_results.data, file = "results/protein_coding/diag/tau_pc_meta_results.Rdata")

# AB
ab_pc_meta_p_values.list <- pbmclapply(genes_adni, function(gene) {
  skat_o_meta(pca_adni.data,"RC2",pca_emif.data,"RC2",csf_cov_diag_adni.data,csf_cov_diag_emif.data,gene)
}, mc.cores = 12)

ab_pc_meta_results.data <- do.call(rbind, ab_pc_meta_p_values.list)
save(ab_pc_meta_results.data, file = "results/protein_coding/diag/ab_pc_meta_results.Rdata")

# NFL
nfl_pc_meta_p_values.list <- pbmclapply(genes_adni, function(gene) {
  skat_o_meta(pca_adni.data,"RC3",pca_emif.data,"RC3",csf_cov_diag_adni.data,csf_cov_diag_emif.data,gene)
}, mc.cores = 12)

nfl_pc_meta_results.data <- do.call(rbind, nfl_pc_meta_p_values.list)
save(nfl_pc_meta_results.data, file = "results/protein_coding/diag/nfl_pc_meta_results.Rdata")

# YKL40
ykl40_pc_meta_p_values.list <- pbmclapply(genes_adni, function(gene) {
  skat_o_meta(pca_adni.data,"RC4",pca_emif.data,"RC4",csf_cov_diag_adni.data,csf_cov_diag_emif.data,gene)
}, mc.cores = 12)

ykl40_pc_meta_results.data <- do.call(rbind, ykl40_pc_meta_p_values.list)
save(ykl40_pc_meta_results.data, file = "results/protein_coding/diag/ykl40_pc_meta_results.Rdata")

# Ng
ng_pc_meta_p_values.list <- pbmclapply(genes_adni, function(gene) {
  skat_o_meta(pca_adni.data,"RC5",pca_emif.data,"RC5",csf_cov_diag_adni.data,csf_cov_diag_emif.data,gene)
}, mc.cores = 12)

ng_pc_meta_results.data <- do.call(rbind, ng_pc_meta_p_values.list)
save(ng_pc_meta_results.data, file = "results/protein_coding/diag/ng_pc_meta_results.Rdata")


