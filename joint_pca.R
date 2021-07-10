library(missMDA)
library(psych)
library(dplyr)

# Load INT transformed CSF biomarkers
load("phenotypes/adni/csf_adni.Rdata")
load("phenotypes/emif/csf_emif.Rdata")

# Add study identifier
csf_adni.data$study <- "adni"
csf_emif.data$study <- "emif"

# Reduce to relevant variables
adni_biomarkers <- c("TAU","PTAU","CSFNFL","ABETA","YKL40","CSFNG")
emif_biomarkers <- c("Ttau_ASSAY_Zscore","Ptau_ASSAY_Zscore","Central_CSF_NFL","Central_CSF_AB42","Central_CSF_YKL40","Central_CSF_Neurogranin")

adni_biomarkers.data <- csf_adni.data[c("study","PTID",adni_biomarkers)]
emif_biomarkers.data <- csf_emif.data[c("study","Gentli_ID",emif_biomarkers)]

# Harmonize CSF biomarker names
names(adni_biomarkers.data) <- c("study","ID",adni_biomarkers)
names(emif_biomarkers.data) <- c("study","ID",adni_biomarkers)

# Make direction in EMIF for Tau more conventional
emif_biomarkers.data$TAU <- -1*emif_biomarkers.data$TAU
emif_biomarkers.data$PTAU <- -1*emif_biomarkers.data$PTAU

# Check correlations
cor(adni_biomarkers.data[adni_biomarkers], use = "pairwise.complete.obs")
cor(emif_biomarkers.data[adni_biomarkers], use = "pairwise.complete.obs")

# Merge both studies
adni_emif.data <- as.data.frame(rbind(adni_biomarkers.data, emif_biomarkers.data))
adni_emif.data$study <- as.factor(adni_emif.data$study)

# Descriptives
multi.hist(adni_emif.data[adni_biomarkers])
describe(adni_emif.data)

# Use cross-validation to determine optimum number of componencts
set.seed(20200617)
nb <- estim_ncpPCA(adni_emif.data[adni_biomarkers], method.cv = "loo")
nb$ncp
png(file="figures/pca/elbow_plot.png")
plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")
dev.off()

# Perform PCA based imputation
set.seed(20200617)
adni_emif_int.imputePCA <- imputePCA(adni_emif.data[adni_biomarkers], ncp = 5)
adni_emif_int_imputed.data <- as.data.frame(adni_emif_int.imputePCA$completeObs)
adni_emif_int_imputed.data$study <- adni_emif.data$study
adni_emif_int_imputed.data$ID <- adni_emif.data$ID

# PCA
pca.fit <- pca(adni_emif_int_imputed.data[adni_biomarkers], nfactors = 5)
# Save loadings
write.csv(unclass(pca.fit$loadings), "results/pca_biomarkers_int_adni_emif.csv")
# compute PC scores
pca_adni_emif.data <- data.frame(predict(pca.fit, data = adni_emif_int_imputed.data[adni_biomarkers]))
pca_adni_emif.data$study <- adni_emif.data$study
pca_adni_emif.data$ID <- adni_emif.data$ID
# Check PC distribution
multi.hist(pca_adni_emif.data)

# check data
describe(pca_adni_emif.data)

# Save PC data
save(pca_adni_emif.data, file = "phenotypes/joint/pca_adni_emif.Rdata")

# Split per cohort and rename ID variable
pca_adni.data <- pca_adni_emif.data[pca_adni_emif.data$study == "adni", c("ID","RC1","RC2","RC3","RC4","RC5")]
names(pca_adni.data) <- c("PTID","RC1","RC2","RC3","RC4","RC5")
save(pca_adni.data, file = "phenotypes/adni/pca_adni.Rdata")

# Split per cohort and rename ID variable
pca_emif.data <- pca_adni_emif.data[pca_adni_emif.data$study == "emif", c("ID","RC1","RC2","RC3","RC4","RC5")]
names(pca_emif.data) <- c("Gentli_ID","RC1","RC2","RC3","RC4","RC5")
save(pca_emif.data, file = "phenotypes/emif/pca_emif.Rdata")

### Heatmap
library(corrplot)

loadings.mat <- t(unclass(pca.fit$loadings))
rownames(loadings.mat) <- c("Tau pathology/degeneration","Injury/inflammation","Aβ pathology","Non-AD inflammation","Non-AD synaptic functioning")
colnames(loadings.mat) <- c("Tau","pTau","NFL","Amyloid","YKL-40","Ng")

col3 <- colorRampPalette(c("red", "white", "blue")) 

png(file="figures/pca/heatmap.png", width = 4000, height = 3000, pointsize = 80)
corrplot(loadings.mat, tl.srt = 0, col = col3(20), method = "color")
dev.off()

loadings_cor.mat <- cbind(loadings.mat,c(-0.14,-0.20,0.28,-0.02,0.01))
colnames(loadings_cor.mat) <- c("Tau","pTau","NFL","Amyloid","YKL-40","Ng","MMSE cor")

png(file="figures/pca/heatmap_cor.png", width = 4000, height = 3000, pointsize = 80)
corrplot(loadings_cor.mat, tl.srt = 0, col = col3(20), method = "color")
dev.off()

### Residualize
### ADNI
# Load covariates
load("phenotypes/adni/csf_cov_adni.Rdata")
# Merge covariate data
pca_resids_adni.data <- merge(pca_adni.data, csf_cov_adni.data, by = "PTID")

# Residualize
compute_residuals <- function(mediator_name) {
  scale(resid(lm(paste0(mediator_name, " ~ AGE + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"), data = pca_resids_adni.data, na.action=na.exclude)))
}

outcomes <- c(paste0("RC",1:5))
outcomes_resids <- paste0(outcomes, "_resids")
pca_resids_adni.data[outcomes_resids] <- lapply(outcomes, compute_residuals)

save(pca_resids_adni.data, file = "phenotypes/emif/pca_resids_adni.Rdata")

### EMIF
# Load covariates
load("phenotypes/emif/csf_cov_emif.Rdata")
# Merge covariate data
pca_resids_emif.data <- merge(pca_emif.data, csf_cov_emif.data, by = "Gentli_ID")

# Residualize
compute_residuals <- function(mediator_name) {
  scale(resid(lm(paste0(mediator_name, " ~ ageatcsfpet + gender + PC1 + PC2 + PC3 + PC4"), data = pca_resids_emif.data, na.action=na.exclude)))
}

outcomes <- c(paste0("RC",1:5))
outcomes_resids <- paste0(outcomes, "_resids")
pca_resids_emif.data[outcomes_resids] <- lapply(outcomes, compute_residuals)

save(pca_resids_emif.data, file = "phenotypes/emif/pca_resids_emif.Rdata")

### Correlations between PCA and MMSE
# Load and merge residualized MMSE data
load("phenotypes/adni/mmse_resids_adni.Rdata")
pca_mmse_resids_adni.data <- merge(pca_resids_adni.data, adni_resids.data[c("PTID","mmse_residuals")], by = "PTID", all.x = T)

load("phenotypes/emif/mmse_resids_emif.Rdata")
pca_mmse_resids_emif.data <- merge(pca_resids_emif.data, emif_resids.data[c("PTID","mmse_residuals")], by.x = "Gentli_ID", by.y = "PTID", all.x = T)

adni_emif_resids.data <- bind_rows(pca_mmse_resids_adni.data[c("mmse_residuals", outcomes_resids)], pca_mmse_resids_emif.data[c("mmse_residuals",outcomes_resids)])
adni_emif_resids.data$ID <- pca_adni_emif.data$ID
save(adni_emif_resids.data, file = "phenotypes/joint/adni_emif_resids.Rdata")

# Correlation ADNI
round(cor(adni_emif_resids.data[1:127, c("mmse_residuals", outcomes_resids)], use = "complete.obs", method = "spearman"), 2)

# Correlation EMIF
round(cor(adni_emif_resids.data[128:480, c("mmse_residuals", outcomes_resids)], use = "complete.obs", method = "spearman"), 2)

# Correlation Meta
round(cor(adni_emif_resids.data[c("mmse_residuals", outcomes_resids)], use = "complete.obs", method = "spearman"), 2)
summary(lm(mmse_residuals ~ RC3_resids, data = adni_emif_resids.data, na.action=na.exclude))

# PCA ADNI
pca.fit <- pca(adni_emif_int_imputed.data[adni_emif_int_imputed.data$study == "adni",adni_biomarkers], nfactors = 5)
# Save loadings
write.csv(unclass(pca.fit$loadings), "results/pca_biomarkers_int_adni.csv")

# PCA EMIF
pca.fit <- pca(adni_emif_int_imputed.data[adni_emif_int_imputed.data$study == "emif",adni_biomarkers], nfactors = 5)
# Save loadings
write.csv(unclass(pca.fit$loadings), "results/pca_biomarkers_int_emif.csv")

