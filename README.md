###
### Documents related to the assessment of side-effects from batch effect adjustments using
### ComBat with study group as covariate or similar methods.
###

README for combat warning paper.
----------------------

xxxx
### reanalysis/  

The individual public datasets that are used to illustrate the adverse effect of ComBat.
Consist of a figurescript that loads the data, process with or without ComBat, run a significance analysis and creates a plot of the p-values.
For the specially interested there is a much longer reproduction-report in html made from r-markdown (need to be run individually) describing the choices done, limitations and some additional sanity checks.


### commonscripts

- helperfunctions.r  methods not related to a individual data sets
- introductionplot.r Creates figure 1 in the article. Generates the data and doing the different adjustments and plots the results.