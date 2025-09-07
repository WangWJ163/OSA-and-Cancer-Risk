library(TwoSampleMR)
library(vroom)

###1.Read exposure data
exp<-read.csv("Obs.csv")

###2.Read outcome data

out<-vroom('breast.meta.gz',col_names = TRUE)
colnames(out)##Check column names
out = out[,c(4,5,8,10,13)]
colnames(out)
colnames(out) = c("effect_allele","other_allele","p","OR","SNP")
out$beta<-log(out$OR)
out$se<-abs(log(out$OR)/qnorm(out$p/2))
save(out,file="breast cancer.rda")
load("breast cancer.rda")
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
data<-read.csv("harmonise_data.csv")

###6 MR analyses
methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_raps")
results <- mr(data, method_list = methods)
write.csv(results,'results.csv')
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

###10. leave one out analyses

# Get IVW method results
full_model <- mr(data, method_list = "mr_ivw")
ivw_result <- full_model[full_model$method == "Inverse variance weighted", ]

# Extract b and se
b <- ivw_result$b
se <- ivw_result$se

# Convert to OR and 95% CI
or <- exp(b)
or_lower <- exp(b - 1.96 * se)
or_upper <- exp(b + 1.96 * se)

# Create results dataframe
or_results <- data.frame(
  Method = "IVW",
  OR = round(or, 3),
  CI_Lower = round(or_lower, 3),
  CI_Upper = round(or_upper, 3),
  P_Value = round(ivw_result$pval, 4),
  stringsAsFactors = FALSE
)

# Display results
cat("IVW Results as Odds Ratios:\n")
cat("OR:", or_results$OR, "\n")
cat("95% CI: [", or_results$CI_Lower, ", ", or_results$CI_Upper, "]\n")
cat("P-value:", or_results$P_Value, "\n")

# If needed, also convert all leave-one-out analysis results to OR
leave_one_out <- mr_leaveoneout(data)

# Convert each result of the leave-one-out analysis
loo_or <- exp(leave_one_out$b)
loo_or_lower <- exp(leave_one_out$b - 1.96 * leave_one_out$se)
loo_or_upper <- exp(leave_one_out$b + 1.96 * leave_one_out$se)

# Create leave-one-out OR results dataframe with P values
loo_or_results <- data.frame(
  SNP = leave_one_out$SNP,
  OR = round(loo_or, 3),
  CI_Lower = round(loo_or_lower, 3),
  CI_Upper = round(loo_or_upper, 3),
  P_Value = round(leave_one_out$p, 4),  # 添加P值
  stringsAsFactors = FALSE
)

# Add full model OR result with P value
full_or_row <- data.frame(
  SNP = "Full Model (IVW)",
  OR = round(or, 3),
  CI_Lower = round(or_lower, 3),
  CI_Upper = round(or_upper, 3),
  P_Value = round(ivw_result$pval, 4)  # 添加P值
)

# Combine all results
all_or_results <- rbind(loo_or_results, full_or_row)

# Add a column to identify if it's the full model
all_or_results$Analysis_Type <- ifelse(all_or_results$SNP == "Full Model (IVW)", 
                                       "Full Model", "Leave-One-Out")

# earrange the column order so that P_Value is last
all_or_results <- all_or_results[, c("SNP", "OR", "CI_Lower", "CI_Upper", "P_Value", "Analysis_Type")]

# save results
write.csv(all_or_results, "leave_one_out_results.csv", row.names = FALSE)
all_or_results<-read.csv("leave_one_out_results.csv")
# Print table
print(all_or_results)

library(ggplot2)
library(dplyr)
###Prepare data for plotting
plot_data <- all_or_results %>%
  mutate(
    SNP = factor(SNP, levels = c(SNP[SNP != "Full Model (IVW)"], "Full Model (IVW)")),
    is_full = ifelse(SNP == "Full Model (IVW)", "Full Model", "Leave-One-Out"),
    # Create annotation text - OR and CI
    or_ci_text = paste0(round(OR, 3), " (", round(CI_Lower, 3), "-", round(CI_Upper, 3), ")"),
    # Create P-value annotation text
    p_text = ifelse(P_Value < 0.001, "<0.001",
                    ifelse(P_Value < 0.01, sprintf("%.3f", P_Value),
                           sprintf("%.3f", P_Value)))
  )

###Ensure correct SNP order and reverse
snp_levels <- rev(levels(plot_data$SNP)) # Reverse order
plot_data$SNP <- factor(plot_data$SNP, levels = snp_levels)

###Get y-axis range (excluding column headers)
y_min <- 0.5 # Start from first data point
y_max <- length(levels(plot_data$SNP)) + 0.5 # End at last data point

###Create forest plot (OR version) with increased row spacing
or_forest_plot <- ggplot(plot_data, aes(x = OR, y = SNP, color = is_full)) +
  
###Add shortened red reference line (x=1.0)
geom_segment(
  aes(x = 1, xend = 1, y = y_min, yend = y_max),
  color = "red", linetype = "dashed", size = 0.8
) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.1, size = 0.8) + # Reduce height to increase row spacing
  geom_point(size = 3) +
  
###Add OR and CI text annotation (right side)
geom_text(aes(x = max(CI_Upper, na.rm = TRUE) * 1.1, label = or_ci_text),
          hjust = 0, size = 4, color = "black") + # Increase text size
  
###Add P-value text annotation (left side)
geom_text(aes(x = min(CI_Lower, na.rm = TRUE) * 0.7, label = p_text),
          hjust = 0.5, size = 4, color = "black") + # Increase text size
  
###Add column headers
annotate("text", x = min(plot_data$CI_Lower, na.rm = TRUE) * 0.7,
         y = length(levels(plot_data$SNP)) + 0.5, # Reduce y-value to bring header closer to first row
         label = "P Value", fontface = "bold", size = 5, hjust = 0.5) + # Increase text size
  
  annotate("text", x = max(plot_data$CI_Upper, na.rm = TRUE) * 1.1,
           y = length(levels(plot_data$SNP)) + 0.5, # Reduce y-value to bring header closer to first row
           label = "Odds Ratio", fontface = "bold", size = 5, hjust = 0) + # Increase text size
  
  scale_color_manual(values = c("Full Model" = "blue", "Leave-One-Out" = "black")) +
  scale_x_log10(
    breaks = c(0.5, 0.7, 1.0, 1.5, 2.0),
    limits = c(min(plot_data$CI_Lower, na.rm = TRUE) * 0.5,
               max(plot_data$CI_Upper, na.rm = TRUE) * 1.3)
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.08))) + # Reduce top expansion to bring first row closer to title
  labs(
    x = "Odds Ratio",
    y = "Analysis",
    title = "Leave-One-Out Sensitivity Analysis in breast cancer using Pan-cancer data",
    color = "Analysis Type",
    caption = "Error bars represent 95% confidence intervals"
  ) +
  theme_minimal(base_size = 16) + # Increase base font size
  theme(
    text = element_text(family = "Times New Roman"), # Set global font
    legend.position = "bottom",
    legend.text = element_text(size = 16), # Increase legend text size
    legend.title = element_text(size = 16, face = "bold"), # Increase legend title size
    axis.title.x = element_text(size = 20, face = "bold"), # Increase x-axis title size
    axis.title.y = element_text(size = 20, face = "bold"), # Increase y-axis title size
    axis.text.x = element_text(size = 18), # Increase x-axis text size
    axis.text.y = element_text(size = 18, margin = margin(r = 10)), # Increase y-axis text size
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5,
                              margin = margin(b = 0)), # Significantly reduce bottom margin of title
    plot.caption = element_text(size = 16, face = "italic"), # Increase caption size
    plot.margin = margin(10, 10, 10, 10), # Significantly reduce top margin
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", size = 0.2), # Add light gray horizontal grid lines
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), # Set white background
    plot.background = element_rect(fill = "white", color = NA) # Set white background
  )

print(or_forest_plot)

### Save as PDF with same dimensions as TIFF
ggsave("leave_one_out_analysis.pdf", 
       plot = or_forest_plot,
       device = cairo_pdf,  # High-quality PDF engine
       width = 14,          # Same width as your TIFF (in inches)
       height = 10 + nrow(plot_data) * 0.3,  # Same height calculation
       units = "in",        # Using inches to match TIFF dimensions
       dpi = 600)           # Same resolution





