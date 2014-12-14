###
### Functions used in testcase.r
###


# Write to console
concat=function(...,sep='') paste(format(list(...)),sep=sep,collapse=sep)
print.log=function(...,sep='') cat(...,'\n',sep=sep);
print.vars=function(x=NULL,...) {
  if (!is.null(x)) {
    if (is.character(x)) cat(x,'\n')
    else cat(as.character(substitute(x)),'=',c(x),'\n')
    print.vars(...)
  }
}


# Save plot to file
save.plot.enable=FALSE
save.plot.suffix=NULL
save.plot=function(...,folder='Plots',suffix=save.plot.suffix,sep='_') {
  if (!save.plot.enable) return(invisible())
  filename=paste(paste(c(...,suffix),sep=sep,collapse=sep),'.png',sep='')
  if (!is.null(folder)) filename=file.path(folder,filename)
  print.log('Save plot to ',filename)
  dev.copy(png,filename)
  dev.off()
  return(invisible())
}


# Draw abline and store description for legend
Line.reset=function() {
  Line.leg<<-c();
  Line.col<<-c();
}
Line.add=function(x0,slope,col='gray',lwd=1,label=NULL) {
  if (is.null(label)) label=slope;
  Line.leg<<-c(Line.leg,label);
  Line.col<<-c(Line.col,col);
  abline(x0,slope,col=col,lwd=lwd);
}


# Random matrix and constant matrix
rmatrix=function(rownames,colnames,rfunc=rnorm,...) {
  X=matrix(rfunc(length(rownames)*length(colnames),...),
      length(rownames),length(colnames));
  rownames(X)=rownames;
  colnames(X)=colnames;
  return(X);
}
cmatrix=function(rownames,colnames,value=NA) {
  return(rmatrix(rownames,colnames,rfunc=function(n) rep(value,n)));
}


# RandomData object:
#   N = list of vectors giving 
RandomData=function(N,N.genes,sd.error=1,sd.batch=1,desc=NULL) {
  obj=new.env();
  class(obj)='RandomData';
  obj$N=N;
  obj$N.genes=N.genes;
  obj$sd.error=sd.error;
  obj$sd.batch=sd.batch;
  obj$desc=c(desc,attr(N,'desc'),paste0('err',100*sd.error/sd.batch,'pct'));

  obj$N.samples=sum(sapply(obj$N,sum));
  obj$N.groups=max(sapply(obj$N,length));
  obj$N.batches=length(obj$N);
  obj$samplenames=paste('Sample',1:obj$N.samples);
  obj$batchnames=paste('Batch',1:obj$N.batches);
  obj$groupnames=paste('Group',1:obj$N.groups);
  obj$genenames=paste('Gene',1:obj$N.genes);

  obj$pheno=data.frame(cmatrix(obj$samplenames,c('Batch','Group')));
  obj$pheno[,1]=rep(obj$batchnames,sapply(obj$N,sum));
  obj$pheno[,2]=unlist(lapply(obj$N,function(n) rep(obj$groupnames,n)));

  obj$edata.noise=rmatrix(obj$genenames,obj$samplenames,sd=obj$sd.error);
  obj$effects.batch=rmatrix(obj$genenames,obj$batchnames,sd=obj$sd.batch);
  obj$edata.batch=cmatrix(obj$genenames,obj$samplenames);
  for (sample in 1:obj$N.samples) {
    obj$edata.batch[,sample]=obj$effects.batch[,obj$pheno[sample,'Batch']];
  }
  obj$edata=obj$edata.batch+obj$edata.noise;

  obj$N.batch.effective=sapply(obj$N,effective.N);
  obj$N.effective=sum(obj$N.batch.effective);
  obj$N.batch.size=sapply(obj$N,sum);
  obj$N.group.size=Reduce("+",obj$N);
  obj$N.effective.group=effective.N(obj$N.group.size);
  obj$N.adjust0=obj$N.effective/obj$N.effective.group;
  obj$N.adjust=ComBat.Fadjust(obj,(obj$sd.batch/obj$sd.error)^2);

  obj$ComBat=function(mod) {
    obj$edata.combat=ComBat(dat=obj$edata,batch=obj$pheno$Batch,mod=mod,numCovs=NULL,par.prior=TRUE,prior.plots=FALSE);
    return(obj$edata.combat);
  }

  return(obj);
}


# Compute F adjustment factor for data with shrinked batch adjustments
ComBat.Fadjust=function(obj,nu) {
	obj$nu=nu;
	with(obj,{
		obj$N.AB=matrix(unlist(N),N.groups,N.batches);
		obj$N.A=diag(N.group.size);
		obj$N.B=diag(N.batch.size);
		obj$N.full=rbind(cbind(N.A,N.AB),cbind(t(N.AB),N.B));
		obj$N.A.inv=solve(N.A);
		obj$Dv=diag(1-1/(1+nu*N.batch.size));
		obj$Ev=diag(1/(1+nu*N.batch.size));
		e=matrix(c(rep(1,N.groups),rep(-1,N.batches)));
		E=e%*%t(e);
		obj$U=Ev%*%t(N.AB)%*%N.A.inv;
		IU=rbind(diag(N.groups),U);
		obj$XX.inv=ginv(N.full);
		XX.inv.B=ginv(N.B-t(N.AB)%*%N.A.inv%*%N.AB);
		obj$Var.a=t(IU)%*%XX.inv%*%IU;
		obj$adjust=alttr(N.A.inv)/alttr(Var.a);
		### Alternative computation
		#V=N.A.inv %*% N.AB %*% Dv;
		#Var.a.alt=N.A.inv + V %*% X.inv.B %*% t(V);
		#adjust=alttr(N.A.inv)/alttr(Var.a.alt);
	});
	return(obj$adjust);
}


# Run F test and store results in an object
Ftest = function (dat, mod, mod0,adjust=1)
{
  obj=new.env()
  class(obj)="Ftest"
  obj$arg.adjust=adjust
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0, m)
  Id <- diag(n)
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  rss1 <- rowSums(resid * resid)
  rm(resid)
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*%  t(mod0))
  rss0 <- rowSums(resid0 * resid0)
  rm(resid0)
  fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
  if (is.null(adjust)) {
    adjust=mean(rss1/(n - df1))/mean((rss0 - rss1)/(df1 - df0))
  }
  fadj=fstats*adjust
  p <- 1 - pf(fadj, df1 = (df1 - df0), df2 = (n - df1))
  obj$n1=n; obj$n0=m; obj$df1=df1-df0; obj$df0=n-df1;
  obj$F0=fstats; obj$F=fadj;
  obj$SS0=rss1; obj$MS0=rss1/(n-df1);
  obj$SS1=rss0-rss1; obj$MS1=(rss0-rss1)/(df1-df0);
  obj$P=p;
  obj$adjust=adjust;
  obj$hist=function(m=50) {
    if (is.null(obj$arg.adjust)) xlabel=concat('P values from empirically adjusted ANOVA (adjust=',obj$adjust,')',sep='')
    else if (obj$arg.adjust==1) xlabel='P values from ANOVA'
    else xlabel=concat('P values from adjusted ANOVA (adjust=',obj$adjust,')',sep='');
    hist(obj$P,breaks=(0:m)/m,xlab=xlabel,main=paste('Histogram of P-values'))
  };
  return(obj);
}


# Set up multi*m+full batches with m groups
# - for each group, add multi batches with n0 samples in one group and n1 in the rest
# - add full batches with n0 samples from each group
UnbalancedBatches=function(m,n0,n1=0,full=0,multi=1) {
  N=rep(list(rep(n1,m)),m);
  for (i in 1:m) N[[i]][i]=n0;
  N=c(rep(list(rep(n0,m)),full),rep(N,multi));
  desc=NULL
  if (multi>0) desc=c(desc,paste0(multi,'x',m,'x',n0,'+',n1))
  if (full>0) desc=c(desc,paste0(full,'x',n0))
  attr(N,'desc')=desc
  return(N);
}


# Matrix trace
tr=function(m) sum(diag(m));
alttr=function(m) tr(m)-mean(m)*dim(m)[1];

# Matrix projections
projection=function(m) m%*%ginv(t(m)%*%m)%*%t(m);
perpendicular=function(m) diag(nrow(m))-projection(m);

# Compute effective number of samples in the sense of
#   n_eff = n * sum_i p_i*(1-p_i)
effective.N=function(x) sum(x)-sum(x*x)/sum(x); 

