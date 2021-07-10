library(data.table)
library(qqman)
library(tidyr)

# Function to load data to a specific name from https://stackoverflow.com/a/25455968
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load effective snp number/carriers
load("results/protein_coding/nodiag/protein_coding_annot.Rdata")
load("results/lof/nodiag/lof_annot.Rdata")

# Load annotation
annotation <- fread("../adni_emif_genotypes/annotations/nopatches_biomartresults.tsv")
annotation <- annotation[ ,.(SetID, Chromosome, GeneStart)]

# Define function to summarize results
post_meta_protein <- function(summary, out, out_qq, out_manhattan) {
  # Load summary statistics
  summary.data <- loadRData(summary)
  
  # Define gene names as character (instead of factor)
  summary.data$gene_name_adni <- as.character(summary.data$gene_name_adni)
  
  # Replace snp and carrier counts with effective counts
  summary.data[,c("snp_count_adni","n_carriers_adni","snp_count_emif","n_carriers_emif","snp_total","snp_overlap")] <- protein_coding_annot.data[,c("snp_count_adni","n_carriers_adni","snp_count_emif","n_carriers_emif","snp_total","snp_overlap")]
  
  # Remove genes with missing data, implausibly low p-values and at least 1 carrier in adni
  summary.data <- summary.data[!is.na(summary.data$p.value_adni) &
                                 summary.data$p.value_adni > 1E-100 & 
                                 summary.data$n_carriers_adni > 1, ]
  
  # Remove empty rows
  summary.data <- drop_na(summary.data)
  
  # Merge with annoutation data
  summary.data <- merge(summary.data, annotation, by.x = "gene_name_adni", by.y = "SetID")
  
  # Calculate genome-wide significant and suggestive threshold
  n_genes <- dim(summary.data)[1]
  bonferroni <- 0.05/n_genes
  suggestive <- 1/n_genes
  
  # Write genome-wide significant results
  bonferroni.data <- summary.data[summary.data$p.value_adni < bonferroni, ]
  write.csv(bonferroni.data, file = out, row.names = F, quote = F)
  
  # Calculate lambda
  z=qnorm(summary.data$p.value_adni/2)
  lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
  
  # QQ-plot
  png(file=out_qq, width = 4000, height = 4000, pointsize = 100)
  qq(summary.data$p.value_adni)
  text(x = 1, y = 3, paste("λ =", round(lambda, 2)))
  dev.off()
  
  # Manhattan plot
  # If hits present
  if(dim(bonferroni.data)[1] > 0) {
    png(file=out_manhattan, width = 4000, height = 3000, pointsize = 90)
    manhattan(summary.data, chr = "Chromosome", bp = "GeneStart", p = "p.value_adni", snp = "gene_name_adni", genomewideline = -log10(bonferroni), suggestiveline = -log10(suggestive), annotatePval = bonferroni)
    dev.off()
  } else {
    # If no hits
    png(file=out_manhattan, width = 4000, height = 3000, pointsize = 90)
    manhattan(summary.data, chr = "Chromosome", bp = "GeneStart", p = "p.value_adni", genomewideline = -log10(bonferroni), suggestiveline = -log10(suggestive))
    dev.off()
  }
}

# Protein-coding no diagnosis adjustment
post_meta_protein("results/protein_coding/nodiag/ab_pc_meta_results.Rdata","results/protein_coding/nodiag/ab_pc_hits_adni.csv","figures/protein_coding/nodiag/qq_ab_pc_meta_results_adni.png","figures/protein_coding/nodiag/manhattan_ab_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/nodiag/tau_pc_meta_results.Rdata","results/protein_coding/nodiag/tau_pc_hits_adni.csv","figures/protein_coding/nodiag/qq_tau_pc_meta_results_adni.png","figures/protein_coding/nodiag/manhattan_tau_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/nodiag/nfl_pc_meta_results.Rdata","results/protein_coding/nodiag/nfl_pc_hits_adni.csv","figures/protein_coding/nodiag/qq_nfl_pc_meta_results_adni.png","figures/protein_coding/nodiag/manhattan_nfl_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/nodiag/ykl40_pc_meta_results.Rdata","results/protein_coding/nodiag/ykl40_pc_hits_adni.csv","figures/protein_coding/nodiag/qq_ykl40_pc_meta_results_adni.png","figures/protein_coding/nodiag/manhattan_ykl40_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/nodiag/ng_pc_meta_results.Rdata","results/protein_coding/nodiag/ng_pc_hits_adni.csv","figures/protein_coding/nodiag/qq_ng_pc_meta_results_adni.png","figures/protein_coding/nodiag/manhattan_ng_pc_meta_results_adni.png")

# Protein-coding diagnosis adjustment
post_meta_protein("results/protein_coding/diag/ab_pc_meta_results.Rdata","results/protein_coding/diag/ab_pc_hits_adni.csv","figures/protein_coding/diag/qq_ab_pc_meta_results_adni.png","figures/protein_coding/diag/manhattan_ab_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/diag/tau_pc_meta_results.Rdata","results/protein_coding/diag/tau_pc_hits_adni.csv","figures/protein_coding/diag/qq_tau_pc_meta_results_adni.png","figures/protein_coding/diag/manhattan_tau_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/diag/nfl_pc_meta_results.Rdata","results/protein_coding/diag/nfl_pc_hits_adni.csv","figures/protein_coding/diag/qq_nfl_pc_meta_results_adni.png","figures/protein_coding/diag/manhattan_nfl_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/diag/ykl40_pc_meta_results.Rdata","results/protein_coding/diag/ykl40_pc_hits_adni.csv","figures/protein_coding/diag/qq_ykl40_pc_meta_results_adni.png","figures/protein_coding/diag/manhattan_ykl40_pc_meta_results_adni.png")
post_meta_protein("results/protein_coding/diag/ng_pc_meta_results.Rdata","results/protein_coding/diag/ng_pc_hits_adni.csv","figures/protein_coding/diag/qq_ng_pc_meta_results_adni.png","figures/protein_coding/diag/manhattan_ng_pc_meta_results_adni.png")

# Define function to summarize results
post_meta_lof <- function(summary, out, out_qq, out_manhattan) {
  # Load summary statistics
  summary.data <- loadRData(summary)
  
  # Define gene names as character (instead of factor)
  summary.data$gene_name_adni <- as.character(summary.data$gene_name_adni)
  
  # Replace snp and carrier counts with effective counts
  summary.data[,c("snp_count_adni","n_carriers_adni","snp_count_emif","n_carriers_emif","snp_total","snp_overlap")] <- lof_annot.data[,c("snp_count_adni","n_carriers_adni","snp_count_emif","n_carriers_emif","snp_total","snp_overlap")]
  
  # Remove genes with missing data, implausibly low p-values and at least 1 carrier in adni
  summary.data <- summary.data[!is.na(summary.data$p.value_adni) &
                                 summary.data$p.value_adni > 1E-100 & 
                                 summary.data$n_carriers_adni > 1, ]
  
  # Remove empty rows
  summary.data <- drop_na(summary.data)
  
  # Merge with annoutation data
  summary.data <- merge(summary.data, annotation, by.x = "gene_name_adni", by.y = "SetID")
  
  # Calculate genome-wide significant and suggestive threshold
  n_genes <- dim(summary.data)[1]
  bonferroni <- 0.05/n_genes
  suggestive <- 1/n_genes
  
  # Write genome-wide significant results
  bonferroni.data <- summary.data[summary.data$p.value_adni < bonferroni, ]
  write.csv(bonferroni.data, file = out, row.names = F, quote = F)
  
  # Calculate lambda
  z=qnorm(summary.data$p.value_adni/2)
  lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
  
  # QQ-plot
  png(file=out_qq, width = 4000, height = 4000, pointsize = 100)
  qq(summary.data$p.value_adni)
  text(x = 1, y = 3, paste("λ =", round(lambda, 2)))
  dev.off()
  
  # Manhattan plot
  # If hits present
  if(dim(bonferroni.data)[1] > 0) {
    png(file=out_manhattan, width = 4000, height = 3000, pointsize = 90)
    manhattan(summary.data, chr = "Chromosome", bp = "GeneStart", p = "p.value_adni", snp = "gene_name_adni", genomewideline = -log10(bonferroni), suggestiveline = -log10(suggestive), annotatePval = bonferroni)
    dev.off()
  } else {
    # If no hits
    png(file=out_manhattan, width = 4000, height = 3000, pointsize = 90)
    manhattan(summary.data, chr = "Chromosome", bp = "GeneStart", p = "p.value_adni", genomewideline = -log10(bonferroni), suggestiveline = -log10(suggestive))
    dev.off()
  }
}


# LoF no diagnosis adjustment
post_meta_lof("results/lof/nodiag/ab_pc_meta_results.Rdata","results/lof/nodiag/ab_pc_hits_adni.csv","figures/lof/nodiag/qq_ab_pc_meta_results_adni.png","figures/lof/nodiag/manhattan_ab_pc_meta_results_adni.png")
post_meta_lof("results/lof/nodiag/tau_pc_meta_results.Rdata","results/lof/nodiag/tau_pc_hits_adni.csv","figures/lof/nodiag/qq_tau_pc_meta_results_adni.png","figures/lof/nodiag/manhattan_tau_pc_meta_results_adni.png")
post_meta_lof("results/lof/nodiag/nfl_pc_meta_results.Rdata","results/lof/nodiag/nfl_pc_hits_adni.csv","figures/lof/nodiag/qq_nfl_pc_meta_results_adni.png","figures/lof/nodiag/manhattan_nfl_pc_meta_results_adni.png")
post_meta_lof("results/lof/nodiag/ykl40_pc_meta_results.Rdata","results/lof/nodiag/ykl40_pc_hits_adni.csv","figures/lof/nodiag/qq_ykl40_pc_meta_results_adni.png","figures/lof/nodiag/manhattan_ykl40_pc_meta_results_adni.png")
post_meta_lof("results/lof/nodiag/ng_pc_meta_results.Rdata","results/lof/nodiag/ng_pc_hits_adni.csv","figures/lof/nodiag/qq_ng_pc_meta_results_adni.png","figures/lof/nodiag/manhattan_ng_pc_meta_results_adni.png")

# LoF diagnosis adjustment
post_meta_lof("results/lof/diag/ab_pc_meta_results.Rdata","results/lof/diag/ab_pc_hits_adni.csv","figures/lof/diag/qq_ab_pc_meta_results_adni.png","figures/lof/diag/manhattan_ab_pc_meta_results_adni.png")
post_meta_lof("results/lof/diag/tau_pc_meta_results.Rdata","results/lof/diag/tau_pc_hits_adni.csv","figures/lof/diag/qq_tau_pc_meta_results_adni.png","figures/lof/diag/manhattan_tau_pc_meta_results_adni.png")
post_meta_lof("results/lof/diag/nfl_pc_meta_results.Rdata","results/lof/diag/nfl_pc_hits_adni.csv","figures/lof/diag/qq_nfl_pc_meta_results_adni.png","figures/lof/diag/manhattan_nfl_pc_meta_results_adni.png")
post_meta_lof("results/lof/diag/ykl40_pc_meta_results.Rdata","results/lof/diag/ykl40_pc_hits_adni.csv","figures/lof/diag/qq_ykl40_pc_meta_results_adni.png","figures/lof/diag/manhattan_ykl40_pc_meta_results_adni.png")
post_meta_lof("results/lof/diag/ng_pc_meta_results.Rdata","results/lof/diag/ng_pc_hits_adni.csv","figures/lof/diag/qq_ng_pc_meta_results_adni.png","figures/lof/diag/manhattan_ng_pc_meta_results_adni.png")