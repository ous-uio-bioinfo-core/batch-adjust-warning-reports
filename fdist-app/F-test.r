###
### Compute F-distribution for batch-adjusted data with covariates
###

#source('../commonscripts/theoryfunctions.r');
source('theoryfunctions.r');
library(limma);
library(MASS);
library(sva); # Not in use unless ComBat is used

options(digits=5);

fdistcalc = function(input) {
  # Make sample data using provided design
  design = read.table(textConnection(input$design), sep=",", header=TRUE);
  resamples = input$samples;
  set.seed(input$rngseed);
  dt = RandomData(design,resamples,sd.batch=.5,sd.error=1);
  
  # Set up model
  batch = as.factor(dt$pheno$Batch);
  group = as.factor(dt$pheno$Group);
  mod0 = model.matrix(~1,data=dt$pheno); # Null-model (intercept)
  mod = model.matrix(~1+group,data=dt$pheno); # Model (group effect)
  mod1 = mod[,2:ncol(mod)]; # Group effects only
  mod.batch = model.matrix(~batch,data=dt$pheno); # Batch effect model
  
  # Batch adjust data
  dt$adjdata = removeBatchEffect(dt$edata, batch=batch, design=mod);
  
  # Perform ANOVA to obtain F-statistic
  test = Ftest(dt$adjdata,mod,mod0);
  F0=list(df1=test$df1,df0=test$df0,F=test$F)
  
  # Compute degrees of freedom corrected for batch adjustments
  B.A = t(mod1)%*%perpendicular(mod0)%*%mod1;
  B.AC = t(mod1)%*%perpendicular(mod.batch)%*%mod1;
  M = B.A%*%solve(B.AC);
  Meig = eigen(M,only.values=TRUE)$values;
  M1.est = tr(M);
  M2.est = tr(M%*%M);
  scale.est = M2.est/M1.est;
  Fadj = list(df1=M1.est**2/M2.est,df0=F0$df0,scale=scale.est,F=test$F/scale.est);
  
  # Return results
  ret=list(design=design,resamples=resamples,F0=F0,Fadj=Fadj)
  return(ret);
}

fdisttext = function(x) {
  ret="";
  ret=paste(ret,"<p>Comparison of groups ",paste(colnames(x$design),collapse=','),"</p>");
  ret=paste0(ret,"<p>Standard ANOVA assumes F-statistic ~ F(",x$F0$df1,",",x$F0$df0,")</p>");
  ret=paste0(ret,"<p>Corrected assumption: F-statistic / ",
             x$Fadj$scale," ~ F(",x$Fadj$df1,",",x$Fadj$df0,")</p>");
  return(ret);
}

fdistplot = function(x) {
  N = x$resamples;
  x$F0$quant = qf(ppoints(N),x$F0$df1,x$F0$df0);# Original F-distribution
  x$Fadj$quant = qf(ppoints(N),x$Fadj$df1,x$Fadj$df0);# Original F-distribution
  max.x=sqrt(max(x$F0$quant,x$Fadj$quant))
  max.y=sqrt(max(x$F0$F,x$Fadj$F))
  qqplot(sqrt(x$F0$quant),sqrt(x$F0$F),col='red',
         xlim=c(0,max.x),ylim=c(0,max.y),
         xlab='Quantiles for sqrt(F)',ylab='Standard and corrected sqrt(F)');
  points(qqplot(sqrt(x$Fadj$quant),sqrt(x$Fadj$F),plot.it=FALSE),col='blue');
  abline(0,1,col='gray',lwd=2);
  legend('topleft',cex=1,pch=1,legend=c('Unadjusted (ANOVA)','Corrected'),col=c('red','blue'))
}
