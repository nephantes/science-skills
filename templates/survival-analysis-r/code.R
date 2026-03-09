library(survival)
library(ggplot2)

df <- read.csv("clinical_data.csv")
# Expect columns: time, event (1=event, 0=censored), group

df$group <- factor(df$group)

# Kaplan-Meier fit
km_fit <- survfit(Surv(time, event) ~ group, data = df)
cat("Kaplan-Meier summary:\n")
print(km_fit)

# Extract KM data for ggplot
km_data <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  upper = km_fit$upper,
  lower = km_fit$lower,
  group = rep(levels(df$group), km_fit$strata)
)

p <- ggplot(km_data, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.15, stat = "identity") +
  labs(x = "Time", y = "Survival Probability", title = "Kaplan-Meier Survival Curves") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1))
ggsave("km_curves.png", p, width = 10, height = 6, dpi = 100)

# Log-rank test
lr <- survdiff(Surv(time, event) ~ group, data = df)
cat("\nLog-rank test chi-squared:", round(lr$chisq, 3),
    "p-value:", format.pval(1 - pchisq(lr$chisq, length(lr$n) - 1)), "\n")

# Cox proportional hazards
cox_fit <- coxph(Surv(time, event) ~ group, data = df)
cat("\nCox PH model summary:\n")
print(summary(cox_fit))