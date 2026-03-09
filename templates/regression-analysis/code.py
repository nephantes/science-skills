import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel

df = pd.read_csv("data.csv")
print(f"Data: {df.shape[0]} rows, {df.shape[1]} columns")

# For ordinal logistic regression
y = df['outcome_ordinal']  # ordinal outcome (e.g., severity 1-5)
X = df[['treatment', 'age', 'sex']]
X = pd.get_dummies(X, drop_first=True)

model = OrderedModel(y, X, distr='logit')
result = model.fit(method='bfgs', disp=False)
print(result.summary())
print(f"\nAIC: {result.aic:.2f}")
print(f"BIC: {result.bic:.2f}")

# Odds ratios
odds_ratios = np.exp(result.params)
print(f"\nOdds ratios:\n{odds_ratios}")

# For regular logistic regression
y_binary = (df['outcome'] > 0).astype(int)
X_const = sm.add_constant(X)
logit_model = sm.Logit(y_binary, X_const)
logit_result = logit_model.fit(disp=False)
print(f"\nLogistic regression AIC: {logit_result.aic:.2f}")