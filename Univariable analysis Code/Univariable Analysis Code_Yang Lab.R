library(TwoSampleMR)
library(vroom)
###1.Read exposure data
exp<-read.csv("Obs.csv")

###2.Read outcome data
out<-vroom('all cancer.gz',col_names = TRUE)
##Standardize column names："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"p"，"N"
colnames(out)##Check column names
out = out[,c(2,4,5,7,11,12,13)]
colnames(out)
colnames(out) = c("SNP","effect_allele","other_allele","eaf","beta","se","p")
save(out,file="all cancer.rda")
load("all cancer.rda")

##3.Identify the overlap between exposure and outcome to determine which SNPs are shared by both.
snps<-intersect(exp$SNP,out$SNP)
b<-as.data.frame(cbind(snps,1))
##4.Extract the information of the overlapping SNPs from the exposure and outcome datasets.

exp1<-subset(exp,SNP %in% b$snps)
out1<-subset(out,SNP %in% b$snps)

out1<-format_data(
  out1,
  type = "outcome",
  snps = snps,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
)


###5.Harmonize the effect alleles between the exposure (exp) and outcome (out) datasets.
data<-harmonise_data(exposure_dat = exp1,outcome_dat = out1,action = 3)
data<-subset(data,data$palindromic=="FALSE")
write.csv(data,"harmonise_data.csv")

###6 MR analyses
methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_raps")
results <- mr(data, method_list = methods)
results <-generate_odds_ratios(results)   #Convert to Odds Ratios (OR).  
write.csv(results,'results-OR.csv')

###7.Heterogeneity analysis
heterogeneity <- mr_heterogeneity(data) 
heterogeneity
write.csv(heterogeneity,'heterogeneity.csv')

###8.Horizontal pleiotropy analysis
pleio <- mr_pleiotropy_test(data)
pleio
write.csv(pleio,'pleio.csv')

###9.MRPRESSO analysis
library(MRPRESSO)   
mr_presso(BetaOutcome = 'beta.outcome',
          BetaExposure = 'beta.exposure', 
          SdOutcome = 'se.outcome', 
          SdExposure = 'se.exposure', 
          data = data, OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 5000, seed = NULL)


