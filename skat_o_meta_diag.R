# ### For testing purposes
# phenotype_adni <- pca_adni.data
# phenotype_name_adni <- "RC1"
# phenotype_emif <- pca_emif.data
# phenotype_name_emif <-  "RC1"
# covariates_adni <- csf_cov_diag_adni.data
# covariates_emif <- csf_cov_diag_emif.data
# gene <- genes_adni[[1]]

# Function for SKAT-O meta-analysis
skat_o_meta <- function(phenotype_adni,phenotype_name_adni,phenotype_emif,phenotype_name_emif,covariates_adni,covariates_emif,gene) {
  tryCatch(
    {
      ### ADNI
      # Load gene
      genotype_adni <- fread(paste0("../adni_emif_genotypes/adni/protein_coding/plink/",gene))
      # Get covariate names
      covariate_names_adni <- names(covariates_adni[2:15])
      # Get SNP names
      snp_names_adni <- names(genotype_adni)
      `%notin%` <- Negate(`%in%`)
      snp_names_adni <- snp_names_adni[snp_names_adni %notin% c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
      
      # Merge phenotype, covariate and snp names
      phenotype_covariate_adni <- merge(phenotype_adni[c("PTID",phenotype_name_adni)], covariates_adni, by = "PTID")
      phenotype_covariate_adni <- phenotype_covariate_adni[complete.cases(phenotype_covariate_adni), ]
      phenotype_covariate_genotype_adni <- merge(phenotype_covariate_adni, genotype_adni, by.x = "PTID", by.y = "IID")
      
      # Phenotype vector
      y1 <- as.matrix(phenotype_covariate_genotype_adni[,phenotype_name_adni])
      # Covariate matrix
      x1 <- phenotype_covariate_genotype_adni[,covariate_names_adni]
      x1$Intercept <- 1
      x1 <- as.matrix(x1[c("Intercept", covariate_names_adni)])
      # SNP matrix
      z1 <- as.matrix(phenotype_covariate_genotype_adni[,snp_names_adni,drop=F])
      
      # Number of cases
      # Check first whether effect allele is major or minor
      snp_names_major_adni <- colnames(z1)[colMeans(z1, na.rm = T) >= 1]
      snp_names_minor_adni <- colnames(z1)[colMeans(z1, na.rm = T) < 1]
      # Calculate number of cases, taking into account right coding of effect allele
      if(!identical(snp_names_major_adni, character(0)) & identical(snp_names_minor_adni, character(0))) {
        n_carriers_adni <- sum(rowMeans(z1[,snp_names_major_adni,drop=F], na.rm = T) < 2, na.rm = T)
      }
      if(identical(snp_names_major_adni, character(0)) & !identical(snp_names_minor_adni, character(0))) {
        n_carriers_adni <- sum(rowMeans(z1[,snp_names_minor_adni,drop=F], na.rm = T) > 0, na.rm = T)
      }
      if(!identical(snp_names_major_adni, character(0)) & !identical(snp_names_minor_adni, character(0))) {
        n_carriers_adni <- sum(rowMeans(z1[,snp_names_major_adni,drop=F], na.rm = T) < 2 | rowMeans(z1[,snp_names_minor_adni,drop=F], na.rm = T) > 0, na.rm = T)
      }
      
      # Sample size
      n_adni <- dim(x1)[1]
      
      # Fit null model
      obj.null_adni <-SKAT_Null_Model(y1 ~ x1)
      # Combined SKAT and BURDEN test
      obj_adni = SKAT(z1, obj.null_adni, method = "SKATO", is_dosage = T)
      # Save p-value
      p.value_adni <- obj_adni$p.value
      
      ### EMIF
      # Load gene
      genotype_emif <- fread(paste0("../adni_emif_genotypes/emif/protein_coding/plink/",gene))
      # Get covariate names
      covariate_names_emif <- names(covariates_emif[2:9])
      # Get SNP names
      snp_names_emif <- names(genotype_emif)
      `%notin%` <- Negate(`%in%`)
      snp_names_emif <- snp_names_emif[snp_names_emif %notin% c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
      
      # Merge phenotype, covariate and snp names
      phenotype_covariate_emif <- merge(phenotype_emif[c("Gentli_ID",phenotype_name_emif)], covariates_emif, by = "Gentli_ID")
      phenotype_covariate_emif <- phenotype_covariate_emif[complete.cases(phenotype_covariate_emif), ]
      phenotype_covariate_genotype_emif <- merge(phenotype_covariate_emif, genotype_emif, by.x = "Gentli_ID", by.y = "IID")
      
      # Phenotype vector
      y2 <- as.matrix(phenotype_covariate_genotype_emif[,phenotype_name_emif])
      # Covariate matrix
      x2 <- phenotype_covariate_genotype_emif[,covariate_names_emif]
      x2$Intercept <- 1
      x2 <- as.matrix(x2[c("Intercept", covariate_names_emif)])
      # SNP matrix
      z2 <- as.matrix(phenotype_covariate_genotype_emif[,snp_names_emif,drop=F])
      
      # Number of cases
      # Check first whether effect allele is major or minor
      snp_names_major_emif <- colnames(z2)[colMeans(z2, na.rm = T) >= 1]
      snp_names_minor_emif <- colnames(z2)[colMeans(z2, na.rm = T) < 1]
      # Calculate number of cases, taking into account right coding of effect allele
      if(!identical(snp_names_major_emif, character(0)) & identical(snp_names_minor_emif, character(0))) {
        n_carriers_emif <- sum(rowMeans(z2[,snp_names_major_emif,drop=F], na.rm = T) < 2, na.rm = T)
      }
      if(identical(snp_names_major_emif, character(0)) & !identical(snp_names_minor_emif, character(0))) {
        n_carriers_emif <- sum(rowMeans(z2[,snp_names_minor_emif,drop=F], na.rm = T) > 0, na.rm = T)
      }
      if(!identical(snp_names_major_emif, character(0)) & !identical(snp_names_minor_emif, character(0))) {
        n_carriers_emif <- sum(rowMeans(z2[,snp_names_major_emif,drop=F], na.rm = T) < 2 | rowMeans(z2[,snp_names_minor_emif,drop=F], na.rm = T) > 0, na.rm = T)
      }
      
      # Sample size
      n_emif <- dim(x2)[1]
      
      # Fit null model
      obj.null_emif <-SKAT_Null_Model(y2 ~ x2)
      # Combined SKAT and BURDEN test across all phenotype kernels
      obj_emif = SKAT(z2, obj.null_emif, method = "SKATO", is_dosage = T)
      # Save p-value
      p.value_emif <- obj_emif$p.value
      
      ### Meta-analysis
      # List of outcome per study
      y.list <- list(y1,y2)
      # List of covariate per study
      x.list <- list(x1,x2)
      
      # Merge genotypes from both studyies into one matrix
      z1 <- phenotype_covariate_genotype_adni[,snp_names_adni,drop=F]
      z2 <- phenotype_covariate_genotype_emif[,snp_names_emif,drop=F]
      z <- as.matrix(bind_rows(z1,z2))
      
      # Fit null model
      obj.null <- Meta_Null_Model(y.list, x.list, n.cohort=2, out_type="C")
      # Combined SKAT and BURDEN test across all phenotype kernels
      obj <- MetaSKAT_wZ(z, obj.null, method="optimal", is.separate = TRUE)
      p.value_meta <- obj$p.value
      
      ### Additional information
      # Gene name
      gene_name_adni <- tools::file_path_sans_ext(basename(gene))
      # Sample size
      n_total <- n_adni + n_emif
      # Number of carriers
      n_carriers_total <- n_carriers_adni + n_carriers_emif
      # Number of snps and overlap
      snp_count_adni <- length(snp_names_adni)
      snp_count_emif <- length(snp_names_emif)
      snp_total <- length(unique(c(snp_names_adni,snp_names_emif)))
      snp_overlap <- sum(snp_names_adni %in% snp_names_emif)
      
      # Save results in a data.frame
      result <- data.frame(gene_name_adni,n_adni,snp_count_adni,n_carriers_adni,p.value_adni,n_emif,snp_count_emif,n_carriers_emif,p.value_emif,n_total,snp_total,snp_overlap,n_carriers_total,p.value_meta)
      return(result)
    },
    error=function(e) {
      gene_name_adni <- tools::file_path_sans_ext(basename(gene))
      n_adni <- as.numeric(NA)
      snp_count_adni <- as.numeric(NA)
      n_carriers_adni <- as.numeric(NA)
      p.value_adni <- as.numeric(NA)
      n_emif <- as.numeric(NA)
      snp_count_emif <- as.numeric(NA)
      n_carriers_emif <- as.numeric(NA)
      p.value_emif <- as.numeric(NA)
      n_total <- as.numeric(NA)
      snp_total <- as.numeric(NA)
      snp_overlap <- as.numeric(NA)
      n_carriers_total <- as.numeric(NA)
      p.value_meta <- as.numeric(NA)
      result <- data.frame(gene_name_adni,n_adni,snp_count_adni,n_carriers_adni,p.value_adni,n_emif,snp_count_emif,n_carriers_emif,p.value_emif,n_total,snp_total,snp_overlap,n_carriers_total,p.value_meta)
    }
  )
}

