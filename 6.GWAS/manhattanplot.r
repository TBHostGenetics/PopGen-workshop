# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load GLM association test results (modify filename as needed)
glm_results_file <- "gwasresults.txt"  # Change to your actual file path
df <- read.table(glm_results_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure necessary columns exist: 'CHR' (Chromosome), 'BP' (Base Pair Position), 'P' (P-value)
required_cols <- c("CHR", "BP", "P")
if (!all(required_cols %in% colnames(df))) {
  stop("The input file must contain 'CHR' (Chromosome), 'BP' (Base Pair Position), and 'P' (P-value) columns.")
}

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
