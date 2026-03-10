library(tidyverse)
library(readxl)

# Load variant data
variants <- read_excel("variants.xlsx")

# Filter by quality
filtered <- variants %>%
  filter(FILTER == "PASS")

# Calculate VAF per sample
vaf_summary <- filtered %>%
  group_by(sample, gene) %>%
  summarise(
    mean_vaf = mean(VAF, na.rm = TRUE),
    n_variants = n(),
    .groups = "drop"
  )

# Compare groups
library(rstatix)
stat_result <- vaf_summary %>%
  group_by(gene) %>%
  t_test(mean_vaf ~ group) %>%
  adjust_pvalue(method = "BH")

write.csv(stat_result, "variant_stats.csv", row.names = FALSE)