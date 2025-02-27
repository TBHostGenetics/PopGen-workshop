library(ggplot2)
library(dplyr)

# File paths (edit these if needed)
glm_file <- "gwasresults.PHENO1.glm.logistic.hybrid.adjusted"  # PLINK GWAS output
bim_file <- "PopGen_allSamples_imputed_QC_unrelated_casecontrolonly.bim"  # PLINK BIM file

# Load GWAS results
gwas <- read.table(glm_file, header=TRUE, stringsAsFactors=FALSE)

# Load SNP positions from BIM file
bim <- read.table(bim_file, header=FALSE, stringsAsFactors=FALSE,
                  col.names=c("CHR", "SNP", "CM", "BP", "A1", "A2"))

# Select relevant columns from BIM
bim <- bim %>% select(SNP, CHR, BP)

# Merge GWAS results with BP info
gwas <- merge(gwas, bim, by="SNP")

# Reorder columns for clarity
gwas <- gwas %>% select(SNP, CHR, BP, everything())

# Load GLM association test results
df <- read.table(gwas, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure necessary columns exist: 'CHR' (Chromosome), 'BP' (Base Pair Position), 'P' (P-value)
#required_cols <- c("CHR", "BP", "P")
#if (!all(required_cols %in% colnames(df))) {
  #stop("The input file must contain 'CHR' (Chromosome), 'BP' (Base Pair Position), and 'P' (P-value) columns.")
#}

# Convert chromosome to factor for proper ordering
df$CHR <- as.factor(df$CHR)

# Compute -log10(p-value)
df$logP <- -log10(df$P)

# Define genome-wide significance threshold
threshold <- -log10(5e-8)

# Create chromosome offsets for proper grouping on the x-axis
df <- df %>%
  group_by(CHR) %>%
  mutate(chr_offset = as.numeric(CHR) * 1.2e6,  # Assign a space between chromosomes
         adjusted_BP = BP + chr_offset) %>%
  ungroup()

# Calculate chromosome midpoints for x-axis labels
chrom_midpoints <- df %>%
  group_by(CHR) %>%
  summarize(midpoint = median(adjusted_BP))

# Create Manhattan plot
ggplot(df, aes(x = adjusted_BP, y = logP, color = CHR)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 1) +
  scale_color_manual(values = rep(c("blue", "black"), length(unique(df$CHR)) / 2)) +  # Alternate colors per chromosome
  scale_x_continuous(breaks = chrom_midpoints$midpoint, labels = chrom_midpoints$CHR) +  # Show chromosome labels
  labs(title = "Manhattan Plot of GLM Association Test",
       x = "Chromosome",
       y = expression(-log[10](italic(P)))) +
  theme_minimal() +
  theme(legend.position = "none")  # Hide legend for cleaner visualization
