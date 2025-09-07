library(meta)
# 1. Load data and prepare subset
data <- read.csv("lung cancer_ meta.csv", header = TRUE)
data <- data[c(3, 8), ]  # Ensure this is the correct subset

# 2. Make sure studlab uses data$group (not data$study)
data$study <- data$group  # If group column contains study names

# 3. Perform meta-analysis
metadata <- metagen(
  TE = data$beta, 
  seTE = data$se, 
  studlab = data$study,
  data = data, 
  sm = "OR",
  n.e = data$ncase, 
  n.c = data$ncontrol, 
  pval = data$p,
  random = FALSE, 
  common = TRUE
)

# 4. Generate formatted forest plot with Times New Roman font
pdf("forest_meta_OR-IVW.pdf", width = 12, height = 8, family = "Times")

forest(
  metadata,
  # Title settings
  main = "Meta-Analysis of Lung Cancer Risk",
  main.font = 2,
  main.cex = 1.5,
  print.title = TRUE,
  
  # Column settings - moved OR/CI/P-value to right
  leftcols = c("studlab", "n.e", "n.c", "w.common"),  # Left: Study, Cases, Controls, Weight
  leftlabs = c("Study", "Case", "Control", "Weight"),
  
  rightcols = c("effect", "ci", "pval"),  # Right: OR, 95% CI, P-value
  rightlabs = c("OR", "95% CI", "P value"),
  
  # Graphic styles
  col.square = "blue",
  col.diamond = "red",
  fs.study = 12,
  spacing = 1.2,
  squaresize = 0.5,
  hetstat = TRUE,
  digits = 3,
  digits.pval = 3,
  colgap = "5mm",
  colgap.forest.left = "10mm"
)

dev.off()

