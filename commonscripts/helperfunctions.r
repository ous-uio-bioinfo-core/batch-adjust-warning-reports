

# shuffle samples whithin a batch for selected covariates
# samplenames, vector of strings in order
# batch, vector of strings with the batch id for the samplenames
# covariate, vector of strings with the covariate name for the samplenames
# shufflecovariates, which of the covariate labes that should have their samples shuffeled
# returns the samplenames shuffled.
shufflesamplesinbatch = function(samplenames, batch, covariate, shufflecovariates=NULL)
{
  if(is.null(shufflecovariates))
    shufflecovariates=unique(covariate)
	ret = samplenames
	b = covariate %in% shufflecovariates
	for(thisbatch in unique(batch))
	{
		a = batch== thisbatch		
		ret[a & b] = sample( samplenames[a&b],sum(a&b),  replace=FALSE)		
	}
	return(ret)
}




