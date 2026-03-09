library(multcomp)

df <- read.csv("data.csv")
cat("Data dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Groups:", paste(unique(df$group), collapse=", "), "\n")

# Set control level for Dunnett's test
df$group <- factor(df$group)
df$group <- relevel(df$group, ref = levels(df$group)[1])

# ANOVA
aov_model <- aov(value ~ group, data = df)
cat("\n=== ANOVA ===\n")
print(summary(aov_model))

# Dunnett's test (all vs control)
dunnett <- glht(aov_model, linfct = mcp(group = "Dunnett"))
cat("\n=== Dunnett's Test ===\n")
print(summary(dunnett))
print(confint(dunnett))

# Tukey HSD for all pairwise comparisons
cat("\n=== Tukey HSD ===\n")
print(TukeyHSD(aov_model))