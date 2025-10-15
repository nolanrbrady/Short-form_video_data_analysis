"""
Group Level analysis working on outputs from fnirs_analysis.py
Performs a main effect test then gated post-hoc analysis

Tests:
Using repeated measures ANOVA for main effect
Using FDR for multiple test correction
Using Paired t-test for post-hoc analysis.

Docs referenced:
- Repeated Measures ANOVA: https://www.statsmodels.org/stable/generated/statsmodels.stats.anova.AnovaRM.html#statsmodels-stats-anova-anovarm
- False Detection Rate Correction: https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.fdrcorrection.html#statsmodels.stats.multitest.fdrcorrection
- Paired t-test: https://www.statsmodels.org/stable/generated/statsmodels.stats.weightstats.ttost_paired.html

Author: Nolan Brady
"""