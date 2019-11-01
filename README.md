# The cultural foundations of modern democracies

Data code and instructions to reproduce the findings for the paper "The cultural foundations of modern democracies".

If you make use of any of this code or data please cite: <insert citation>

## directories
data - datasets needed to construct cultural value units 

R - R scripts

python - Python scripts

timeSeriesRegression - time series for each variable in .csv format (used for regressions)


## Get raw data
European Values Survey https://www.gesis.org/en/services/data-analysis/international-survey-programs/european-values-study/

World Values Survey http://www.worldvaluessurvey.org/WVSContents.jsp

## Derive Cosmopolitanism from raw WEVS data 

We provide a file containing the combined World and European Values Survey data called "WEVS"; it contains the 68 common 
cultural value questions sicen 1990, demographic information and the variables are standardized with missing values mean imputed.  Run "ExtractRandC.R" to use weighted principal component analysis to extract Openness to Diversity (C) from WEVS data.

## Derive three Civic confidence variables from WEVS data

Using the file "WEVS_civic", use the script "civicPCA.R" to construct three civic confidence variables: "Institutional Confidence", "Support for Democracy" and "Generalized Trust"

## Show that birth decade differences are independent of time period using model comparison

Run "splitSampesByBirthdecadeAndTimeperiod.py" to split the representative samples for each nation by birth decade and time period for Openness to Diversity and the three Civic Confidence variables. Then "createDataframeForkfold.py" converts all the nation matrices into a dataframe to be used in the model comparison. Then run "kfoldModelComparison.R" which compares hierarchical linear regressions of increasing complexity, testing whether birth decade differences are independent of time period.

## run hierachical time-lagged regression (figure 1)

Compare national time series for Democracy (D), Openness to Diversity  (C) Institutional Confidence (CON), Support for Democracy (S), Generalized Trust (T) and GDP per capita (GDP) (time series provided in the folder "timeSeriesRegression"). Run the file "runRegressions.R" to fit and save results.

Changing the "adultAge" parameter runs regressions assuming an adult age of either 0-10, 10-20 or 20-30 years.  

## run hierachical time-lagged regression (figure 1)

Compare national time series for Democracy (D), Openness to Diversity  (C) Institutional Confidence (CON), Support for Democracy (S), Generalized Trust (T), GDP per capita (GDP) and interaction term UD (time series provided in the folder "timeSeriesRegression"). Run the file "runRegressionsUD.R" to fit and save results.

Changing the "adultAge" parameter runs regressions assuming an adult age of either 0-10, 10-20 or 20-30 years.  

## plot figures 1 and 2 (regression results) 

Run "plotRegressionResults.py" to recreate figure 1 and figure 2.  
