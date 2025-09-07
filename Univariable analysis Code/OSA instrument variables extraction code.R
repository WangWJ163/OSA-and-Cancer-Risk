library(TwoSampleMR)
library(vroom)

## 1. Read exposure data and perform quality control (strong associations + LD pruning) -----------
## ① Read OSA (Obstructive Sleep Apnea) data using vroom
library(vroom)
a <- vroom('finngen_R9_G6_SLEEPAPNO.gz', col_names = TRUE)

## ② Extract required columns and standardize column names ----------
# Standard columns needed: "SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "p", "N"
colnames(a) # Check original column names
a1 = a[,c(3,4,5,7,9,10,11)] # Select relevant columns
colnames(a1) # Verify selected columns
colnames(a1) = c("other_allele","effect_allele","SNP","p","beta","se","eaf") # Rename columns

## ③ Extract genome-wide significant SNPs (p < 5e-8) -------
b <- subset(a1, p < 5e-08) 

## ④ Convert to standardized exposure format ----------
b1 <- format_data(
  b,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "N",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  z_col = "z",
  info_col = "info",
  chr_col = "chr",
  pos_col = "pos",
  log_pval = FALSE
)

## ⑤ Perform LD clumping (online) ----------
exp <- clump_data(b1, clump_kb = 10000, clump_r2 = 0.001)

write.csv(exp, "Obs.csv")

# 2. Calculate F-statistics for instrument strength
dat <- read.csv("Obs.csv")

# Calculate minor allele frequency (MAF)
dat$EAF2 <- (1 - dat$eaf.exposure)
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)

# Function to calculate proportion of variance explained (PVE)
PVEfx <- function(BETA, MAF, SE, N) {
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF)))
  return(pve) 
}

# Apply PVE function to each SNP
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure)

# Calculate F-statistic for each SNP
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))

write.csv(dat, 'Obs.csv')


