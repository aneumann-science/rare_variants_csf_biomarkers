library(SMUT)
library(data.table)
library(dplyr)

load("phenotypes/joint/adni_emif_resids.Rdata")

# Function for mediation analysis protein-coding
mediation_protein <- function(outcome,mediator,gene,data) {
  tryCatch(
    {
      ### ADNI
      # Load gene
      genotype_adni <- fread(paste0("../adni_emif_genotypes/adni/protein_coding/plink/",gene))
      # Get SNP names
      snp_names_adni <- names(genotype_adni)
      `%notin%` <- Negate(`%in%`)
      snp_names_adni <- snp_names_adni[snp_names_adni %notin% c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
      
      ### EMIF
      # Load gene
      genotype_emif <- fread(paste0("../adni_emif_genotypes/emif/protein_coding/plink/",gene))
      # Get SNP names
      snp_names_emif <- names(genotype_emif)
      `%notin%` <- Negate(`%in%`)
      snp_names_emif <- snp_names_emif[snp_names_emif %notin% c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
      
      ### Merge ADNI
      
      # Merge genotypes from both studyies into one matrix
      z1 <- genotype_adni[,..snp_names_adni]
      z2 <- genotype_emif[,..snp_names_emif]
      z <- as.matrix(bind_rows(z1,z2))
      genotype <- as.data.frame(z)
      
      genotype$ID <- c(genotype_adni$IID, genotype_emif$IID)
      
      # Merge genotype, and phenotype
      genotype_phenotype <- merge(genotype,data[c("ID",mediator,outcome)], by = "ID")
      
      # Only keep complete cases of mediator and outcome
      genotype_phenotype <- genotype_phenotype[complete.cases(genotype_phenotype[c(mediator,outcome)]), ]
      
      # Extract genotype matrix and median impute
      z <- as.matrix(genotype_phenotype[c(snp_names_adni,snp_names_emif)[!duplicated(c(snp_names_adni,snp_names_emif))]])
      z_imp <- apply(z, 2, function(x) ifelse(is.na(x), median(x, na.rm=T), x))
      m <- as.matrix(genotype_phenotype[mediator])
      y <- as.matrix(genotype_phenotype[outcome])
      mediation.results <- SMUT(z_imp,m,y)
    },
    error=function(e) {
      print("error")
    }
  )
}

# Function for mediation analysis LoF
mediation_lof <- function(outcome,mediator,gene,data) {
  tryCatch(
    {
      ### ADNI
      # Load gene
      genotype_adni <- fread(paste0("../adni_emif_genotypes/adni/lof/plink/",gene))
      # Get SNP names
      snp_names_adni <- names(genotype_adni)
      `%notin%` <- Negate(`%in%`)
      snp_names_adni <- snp_names_adni[snp_names_adni %notin% c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
      
      ### EMIF
      # Load gene
      genotype_emif <- fread(paste0("../adni_emif_genotypes/emif/lof/plink/",gene))
      # Get SNP names
      snp_names_emif <- names(genotype_emif)
      `%notin%` <- Negate(`%in%`)
      snp_names_emif <- snp_names_emif[snp_names_emif %notin% c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
      
      ### Merge ADNI
      
      # Merge genotypes from both studyies into one matrix
      z1 <- genotype_adni[,..snp_names_adni]
      z2 <- genotype_emif[,..snp_names_emif]
      z <- as.matrix(bind_rows(z1,z2))
      genotype <- as.data.frame(z)
      
      genotype$ID <- c(genotype_adni$IID, genotype_emif$IID)
      
      # Merge genotype, and phenotype
      genotype_phenotype <- merge(genotype,data[c("ID",mediator,outcome)], by = "ID")
      
      # Only keep complete cases of mediator and outcome
      genotype_phenotype <- genotype_phenotype[complete.cases(genotype_phenotype[c(mediator,outcome)]), ]
      
      # Extract genotype matrix and median impute
      z <- as.matrix(genotype_phenotype[c(snp_names_adni,snp_names_emif)[!duplicated(c(snp_names_adni,snp_names_emif))]])
      z_imp <- apply(z, 2, function(x) ifelse(is.na(x), median(x, na.rm=T), x))
      m <- as.matrix(genotype_phenotype[mediator])
      y <- as.matrix(genotype_phenotype[outcome])
      mediation.results <- SMUT(z_imp,m,y)
    },
    error=function(e) {
      print("error")
    }
  )
}

print(mediation_protein("mmse_residuals", "RC3_resids", "IFFO1.raw", adni_emif_resids.data))
print(mediation_protein("mmse_residuals", "RC3_resids", "DTNB.raw", adni_emif_resids.data))
print(mediation_protein("mmse_residuals", "RC3_resids", "NLRC3.raw", adni_emif_resids.data))
print(mediation_protein("mmse_residuals", "RC5_resids", "GABBR2.raw", adni_emif_resids.data))
print(mediation_protein("mmse_residuals", "RC5_resids", "CASZ1.raw", adni_emif_resids.data))
print(mediation_lof("mmse_residuals", "RC3_resids", "SLC22A10.raw", adni_emif_resids.data))
