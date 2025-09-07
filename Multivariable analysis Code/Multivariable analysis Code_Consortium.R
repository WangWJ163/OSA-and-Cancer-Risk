# Multivariable Mendelian Randomization (MVMR) Analysis Pipeline
# Analyzing the effect of multiple exposures on cancer risk

library(tidyr)
library(TwoSampleMR)
library(plyr)
library(MendelianRandomization)
library(MVMR)

# ----------------------------
# 1. DATA LOADING
# ----------------------------

# Load exposure datasets
load("CigarettesPerDay.rda")  # Exposure 1: Cigarettes per day
load("SmokingInitiation.rda") # Exposure 2: Smoking initiation
load("DrinksPerWeek.rda")     # Exposure 3: Alcohol consumption
load("OSA.rda")               # Exposure 4: Obstructive sleep apnea

# ----------------------------
# 2. INSTRUMENT SELECTION
# ----------------------------

## Extract instruments for each exposure (p < 5e-8, LD clumping r2 < 0.001)

# Exposure 1: Cigarettes per day
X1 <- subset(CigarettesPerDay, p < 5e-8)
X1 <- format_data(X1, type = "exposure")
X111 <- clump_data(X1)

# Exposure 2: Smoking initiation
X2 <- subset(SmokingInitiation, p < 5e-8)
X2 <- format_data(X2, type = "exposure")
X222 <- clump_data(X2)

# Exposure 3: Drinks per week
X3 <- subset(DrinksPerWeek, p < 5e-8)
X3 <- format_data(X3, type = "exposure")
X333 <- clump_data(X3)

# Exposure 4: OSA
X4 <- subset(OSA, p < 5e-8)
X4 <- format_data(X4, type = "exposure")
X444 <- clump_data(X4)

# Exposure 5: BMI (extracted from IEU GWAS database)
X555 <- extract_instruments(outcomes = 'ieu-b-4816', p1 = 5e-08, 
                            clump = TRUE, r2 = 0.001, kb = 10000)

# Exposure 6: Diabetes (extracted from EBI GWAS database)
X666 <- extract_instruments(outcomes = 'ebi-a-GCST005413', p1 = 5e-08,
                            clump = TRUE, r2 = 0.001, kb = 10000)

# ----------------------------
# 3. INSTRUMENT PROCESSING
# ----------------------------

# Combine all instruments
SNP_list <- rbind.fill(X111, X222, X333, X444, X555, X666)

# Perform LD clumping on combined instruments
SNP_list_c <- clump_data(SNP_list)

# Extract unique SNPs
SNP <- SNP_list_c$SNP
SNP_u <- unique(SNP)
SNP <- as.data.frame(cbind(SNP_u, 1))
colnames(SNP) <- c("SNP", "nu")
write.csv(SNP, "SNP.csv", row.names = FALSE)

# ----------------------------
# 4. DATA EXTRACTION
# ----------------------------

## Extract SNP information from all exposures and outcome

# Exposure 1: Cigarettes per day
CigarettesPerDay_X1 <- merge(SNP, CigarettesPerDay, by = "SNP")

# Exposure 2: Smoking initiation
SmokingInitiation_X2 <- merge(SNP, SmokingInitiation, by = "SNP")

# Exposure 3: Drinks per week
DrinksPerWeek_X3 <- merge(SNP, DrinksPerWeek, by = "SNP")

# Exposure 4: OSA
OSA_X4 <- merge(SNP, OSA, by = "SNP")

# Exposure 5: BMI
BMI <- extract_outcome_data(snps = SNP$SNP, outcomes = 'ieu-b-4816', proxies = FALSE)
BMI_X5 <- BMI[, c("SNP", "other_allele.outcome", "effect_allele.outcome", 
                  "pval.outcome", "beta.outcome", "se.outcome", "eaf.outcome")]
colnames(BMI_X5) <- c("SNP", "other_allele", "effect_allele", "p", "beta", "se", "eaf")
write.csv(BMI_X5, "BMI_X5.csv")

# Exposure 6: Diabetes
Diabetes <- extract_outcome_data(snps = SNP$SNP, outcomes = 'ebi-a-GCST005413', proxies = FALSE)
Diabetes_X6 <- Diabetes[, c("SNP", "other_allele.outcome", "effect_allele.outcome", 
                            "pval.outcome", "beta.outcome", "se.outcome", "eaf.outcome")]
colnames(Diabetes_X6) <- c("SNP", "other_allele", "effect_allele", "p", "beta", "se", "eaf")
write.csv(Diabetes_X6, "Diabetes_X6.csv")

# Outcome: Cancer
cancer <- extract_outcome_data(snps = SNP$SNP, outcomes = 'ieu-a-1126', proxies = FALSE)
cancer_X7 <- cancer[, c("SNP", "other_allele.outcome", "effect_allele.outcome", 
                        "pval.outcome", "beta.outcome", "se.outcome", "eaf.outcome")]
colnames(cancer_X7) <- c("SNP", "other_allele", "effect_allele", "p", "beta", "se", "eaf")

# ----------------------------
# 5. SNP HARMONIZATION
# ----------------------------

# Find common SNPs across all datasets
snps <- intersect(CigarettesPerDay_X1$SNP, SmokingInitiation_X2$SNP)
snps <- intersect(snps, DrinksPerWeek_X3$SNP)
snps <- intersect(snps, OSA_X4$SNP)
snps <- intersect(snps, BMI_X5$SNP)
snps <- intersect(snps, Diabetes_X6$SNP)
snps <- intersect(snps, cancer_X7$SNP)

# Harmonize effect alleles between exposure and outcome pairs
# (Code continues with harmonization steps as in original script)

# ----------------------------
# 6. MULTIVARIABLE MR ANALYSIS
# ----------------------------

## Using TwoSampleMR package
res <- mv_multiple(mv_DAT)
res_OR <- generate_odds_ratios(res$result)
write.table(res_OR, file = "IVM_OR_TwoSample.xls", sep = "\t", quote = FALSE)

## Using MendelianRandomization package
# Prepare MVMR input
MRMVInputObject <- mr_mvinput(
  bx = cbind(mvmr$betaX1, mvmr$betaX2, mvmr$betaX3, mvmr$betaX4, mvmr$betaX5, mvmr$betaX6),
  bxse = cbind(mvmr$sebetaX1, mvmr$sebetaX2, mvmr$sebetaX3, mvmr$sebetaX4, mvmr$sebetaX5, mvmr$sebetaX6),
  by = mvmr$betaYG,
  byse = mvmr$sebetaYG
)

# IVW method
MRMVObject1 <- mr_mvivw(MRMVInputObject)
mvmr_IVW <- cbind(transform(MRMVObject1@Estimate), transform(MRMVObject1@StdError), transform(MRMVObject1@Pvalue))
colnames(mvmr_IVW) <- c("b", "se", "pval")
OR1 <- generate_odds_ratios(mvmr_IVW)
write.csv(OR1, 'mvmr_IVW_OR.csv', row.names = TRUE)

# Egger method
MRMVObject2 <- mr_mvegger(MRMVInputObject)
mvmr_egger <- cbind(transform(MRMVObject2@Estimate), transform(MRMVObject2@StdError.Est), transform(MRMVObject2@Pvalue.Est))
colnames(mvmr_egger) <- c("b", "se", "pval")
OR2 <- generate_odds_ratios(mvmr_egger)
write.csv(OR2, 'mvmr_egger_OR.csv', row.names = TRUE)

# ----------------------------
# 7. INSTRUMENT STRENGTH & HETEROGENEITY
# ----------------------------

# Using MVMR package
F.data <- format_mvmr(
  BXGs = mvmr[, c(3, 4, 5, 6, 7, 8)],
  BYG = mvmr[, 1],
  seBXGs = mvmr[, c(9, 10, 11, 12, 13, 14)],
  seBYG = mvmr[, 2],
  RSID = "NULL"
)

sres <- strength_mvmr(r_input = F.data, gencov = 0)
write.csv(sres, "F_statistic.csv")



