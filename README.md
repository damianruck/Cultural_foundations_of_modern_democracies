# The cultural foundations of modern democracies

Data code and instructions to reproduce the findings for the paper "The cultural foundations of modern democracies".

If you make use of any of this code or data please cite: *Ruck, D.J., Matthews, L.J., Kyritsis, T. et al. The cultural foundations of modern democracies. Nat Hum Behav (2019).*

## Directories
time_series_normalized - time series for each variable in .csv format (used for regressions)

random_effects - folder contains langauge categories for the multi-level regression

## Get raw data
European Values Survey https://www.gesis.org/en/services/data-analysis/international-survey-programs/european-values-study/

World Values Survey http://www.worldvaluessurvey.org/WVSContents.jsp

## Run hierachical time-lagged regression

Running runRStan.R fits the multilevel time-lagged linear regression to the time series in 'time_series_normalized'.

Compare national time series for Democracy (D), Openness to Diversity  (C) Institutional Confidence (CON), Support for Democracy (S), Generalized Trust (T), GDP per capita (GDP) and interaction term UD (time series provided in the folder "timeSeriesRegression"). Run the file "runRegressionsUD.R" to fit and save results.

Changing the "adultAge" parameter runs regressions assuming an adult age of either 0-10, 10-20 or 20-30 years.  

## Plot figures 1 and 2 (regression results) 

Run "plot.py" to recreate figure 1 and figure 2.  
