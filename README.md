# The cultural foundations of modern democracies

Data code and instructions to reproduce the findings for the paper "The cultural foundations of modern democracies".

If you make use of any of this code or data please cite: <insert citation>

## directories
R - R scripts

python - Python scripts

timeSeriesRegression - time series for each variable in .csv format (used for regressions)

data - additional datasets needed to run results


## Get raw data
European Values Survey https://www.gesis.org/en/services/data-analysis/international-survey-programs/european-values-study/

World Values Survey http://www.worldvaluessurvey.org/WVSContents.jsp


## run hierachical time-lagged regression (figure 1)

Compare national time series for Democracy (D), Openness to Diversity  (C) Institutional Confidence (CON), Support for Democracy (S), Generalized Trust (T) and GDP per capita (GDP) (time series provided in the folder "timeSeriesRegression"). Run the file "runRegressions.R" to fit and save results.

Changing the "adultAge" parameter runs regressions assuming an adult age of either 0-10, 10-20 or 20-30 years.  

## run hierachical time-lagged regression (figure 1)

Compare national time series for Democracy (D), Openness to Diversity  (C) Institutional Confidence (CON), Support for Democracy (S), Generalized Trust (T), GDP per capita (GDP) and interaction term UD (time series provided in the folder "timeSeriesRegression"). Run the file "runRegressionsUD.R" to fit and save results.

Changing the "adultAge" parameter runs regressions assuming an adult age of either 0-10, 10-20 or 20-30 years.  

## plot figures 1 and 2 (regression results) 

Run "plotRegressionResults.py" to recreate figure 1 and figure 2.  
