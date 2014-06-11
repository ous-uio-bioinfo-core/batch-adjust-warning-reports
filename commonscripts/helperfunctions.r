
Xcreatesampleannotation = function( treatment_in_batches)
{
	batches = vector()
	treatments = vector()
	for(i in 1:length(treatment_in_batches))
	{
		#thisbatch = paste("batch", i, sep="")
		thisbatch = i
		for(s in 1:length(treatment_in_batches[[i]]))
		{
			batches = c(batches, rep(thisbatch, treatment_in_batches[[i]][s]))
			#treatments = c(treatments, rep(paste("treatment", s, sep=""), treatment_in_batches[[i]][s]))
			treatments = c(treatments, rep( s, treatment_in_batches[[i]][s]))
		}
	}
	sampleannotation = data.frame(id=paste("sample", 1:length(treatments), sep=""),    treatment = treatments, batch=batches)
	return(sampleannotation)
}




xcreaterandomdata = function(ngenes, sa, mean=0, sd=1)
{
	ret = matrix( rnorm(dim(sa)[1]*ngenes, mean=mean, sd=sd), ncol=dim(sa)[1], nrow=ngenes)
	dimnames(ret) = list( paste("gene", 1:ngenes, sep=""),  sa$id)
	return(ret)
}

xaddbatcheffect = function(df, batches, thismean=1, thissd=1)
{
	ret = df
	for(i in 1:length(unique(batches)))
	{
		for(s in 1:dim(df)[1])
		{
			thisgenesbatcheffect = rnorm(1,mean=thismean, sd=thissd)
			a=batches==i
			#ret[s, a] = ret[s, a] + rnorm(length(ret[s, a]), mean=thisgenesbatcheffect)	
			ret[s, a] = ret[s, a] + thisgenesbatcheffect
		}
	}
	return(ret)
}

# 
xaddconditioneffect = function(df, labels, ndiffgenes, thismean=1, betweensd=1, insidesd=1, affectedconditions)
{
	ret = df
	if(ndiffgenes>0)
	{
		for(i in 1:length(affectedconditions))
		{
			for(s in 1:ndiffgenes)
			{
				thisgeneseffect = rnorm(1, thismean, sd=betweensd)
				a=labels==affectedconditions[i]
				#ret[s, a] = ret[s, a]+ rnorm(length(ret[s, a]), mean=thisgeneseffect)		
				ret[s, a] = rnorm(length(ret[s, a]), mean=thisgeneseffect, sd=insidesd)		
			}
		}
	}
	return(ret)
}


xgetdifftab_sva = function(edata, sa, mod, threshold=0.05)
{
	mod0 = model.matrix(~1,data=sa)
	pValues = f.pvalue(edata,mod,mod0)
	qValues = p.adjust(pValues,method="BH")
	#print(table(qValues<0.05))
	#a = order(pValues)
	#ret = data.frame(gene=names(pValues)[a], p=qValues[a], q=qValues[a])
	ret = data.frame(gene=names(pValues), p=pValues, padjusted=qValues)
	return(ret)
}

#edata = datamatrices[["qnorm"]]
#condition = sampleannotation$covariate1
#contrast = "DP-N"
#block = sampleannotation$chip
#block = NULL

# adapted from
#http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
xgetdifftab_limma = function(edata, condition,  contrast, block=NULL, threshold=0.05)
{
	require(limma)
	fac = factor(condition)	
	design = model.matrix(~0 + fac)
	if(!is.null(block))
	{
		block = as.factor(block)
		#design <- model.matrix(~block+fac)
		design = model.matrix(~0+fac+block)
	}
	colnames(design)=make.names(colnames(design))
	#contrast = paste(unique(sa$treatment)[1], "-", unique(sa$treatment)[2], sep="")
	
	fit <- lmFit(edata, design)
	contrast = paste("fac", contrast, sep="")
	contrast = sub("-", "-fac", contrast)
	cont.matrix = makeContrasts ( contrasts=contrast, levels=design)	
	fit2 = contrasts.fit(fit, cont.matrix)	
	fit2 <- eBayes(fit2)	
	ret = data.frame(gene=names(fit2$Amean), 
                   p=fit2$p.value[,1], 
                   padjusted=p.adjust(fit2$p.value[,1] , method="fdr"),
                   fc=fit2$coefficients[,1],
                   stringsAsFactors=FALSE)
	dimnames(ret)[[1]] = ret$gene 
	return(ret)
}




xpanel.mid <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}



xpanel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}



xgetbalancename = function(sa)
{
	 tab = table(sa[,c("batch", "treatment")])
	 ret = paste("design_" ,  dim(tab)[1], "x", dim(tab)[2], sep="")
	 for(i in 1:dim(tab)[1])
	 {
		thisbatch = paste(tab[i,], collapse="-")
		ret = paste(ret, "_", thisbatch, sep="")
	 }
	 return(ret)
}


#m = combat_mod[i,]
#sa = sampleannotation
#main="after"
xplot_one_gene = function(m, sa, ylim=NA, main="")
{
	a = order(sa$batch)
	
	thispal=c("red", "blue", "green", "orange", "brown", "cyan")
	thissa=sa[a,]
	thisdata=m[a]
	if(any(is.na(ylim)))
	{
		ymax = round(max(thisdata)+2, 1)
		ymin = round(min(thisdata)-2, 1)
		ylim=c(ymin,ymax)
	}
	
	
	plot(thisdata, ylim=ylim, col=thispal[as.factor(thissa$treatment)], main=main,  ylab="Expression")
	
	batchnames=unique(thissa$batch)
	x = match(batchnames, thissa$batch)
	abline(v = x[-1]-0.5)#
	batchlabelx=vector()
	for(i in 1:(length(x)-1))
	{
		batchlabelx = c(batchlabelx, (x[i+1]+x[i])/2)
	}
	batchlabelx = c(batchlabelx, (x[i+1]+length(m))/2)
	#print(x)
	#print(batchlabelx)
	text(labels="Batch", y=ymin+0.5, x=2)
	text(labels=paste(batchnames, sep=""), y=ymin+0.5, x=batchlabelx, cex=1.5)
	#ttestres = t.test(thisdata[sa$])
	treatmentnames=unique(thissa$treatment)
	treatmentmeans = vector()
	for(thistreatment in treatmentnames)
	{
		treatmentmeans[thistreatment]=mean(thisdata[thissa$treatment==thistreatment])
	}
	segments(  x0=(1:length(thisdata))-0.5, y0=treatmentmeans[thissa$treatment], x1=(1:length(thisdata))+0.5, y1=treatmentmeans[thissa$treatment], col=thispal[as.factor(thissa$treatment)])
	legend("topright", legend=paste(treatmentnames, round(treatmentmeans,2)), text.col=thispal[as.factor(treatmentnames)])
}




#dat  = df
#sa = sampleannotation
#name =  paste("sampledata/", getbalancename(sampleannotation), sep="")
xmakecombatinputfiles = function(dat, sa, name)
{
	datafn = paste(name, "_data.txt", sep="")
	write(paste( c("geneinfo", dimnames(dat)[[2]]), collapse="\t"), file=datafn)
	write.table(dat, col.names=FALSE, sep="\t", file=datafn, append=TRUE, quote=FALSE)
	safn= paste(name, "_sampleannotation.txt", sep="")
	ComBatformat = data.frame( sa$id, sa$id, sa$batch, sa$treatment)
	names(ComBatformat)= c("Array name", "Sample name", "Batch", "Covariate 1")
	write.table(ComBatformat, row.names=FALSE, sep="\t", file=safn, quote=FALSE)
}
#makecombatinputfiles(df, sampleannotation, paste("sampledata/", getbalancename(sampleannotation), sep=""))


#samplenames = sampleannotation$code
#batches = sampleannotation$chip
#conditions=sampleannotation$covariate1

# subsetting the samples so that they have a balanced conditon distribution across batches. i.e each batch has the same compitition of conditions (example 70% treated vs 30% untreated)
xdrawbatchbalanceddsamples = function(samplenames, batches, conditions)
{
	# 1. find the minimal count for each condition across batches.
	# 2 find the minimal given count / minimal count inside each batch
	# 3 If greater than 1, scale up all conditions for that batch.
	# 4 draw samples randomly from the above calculated condition-count for each batch
	# 5 create and return the data subset.


	#1
	starttab = table( data.frame(conditions, batches ))	
	tab = starttab[, apply(starttab, 2,  min)>0] # do not use batches with missing conditions
	if(dim(tab)[2]==0)
	{
		warning(  paste("Unable to find a balanced subset that has at least one of each condition. Remove one or more conditions and try again")) 
		return(vector())
	}
	if(dim(tab)[1]<dim(starttab)[1])
	{
		warning(  paste("Lost batches during balancing due to missing conditions in those. ",   "Tried: ", paste( dimnames(starttab)[[2]], collapse=" "), "Got:", paste( dimnames(tab)[[2]], collapse=" "  ) ))
	}
	
	mincond = apply(tab, 1,  min)
	
	#2
	tab2 = mincond/tab
	
	#3
	scalings = apply( 1/tab2, 2, min)
	scaledtab = tab
	scaledtab[,] = mincond	
	scaledtab = t(t(scaledtab) * scalings)
	scaledtab = round(scaledtab,0)
	
	winnersamples= vector()
	for(c in dimnames(scaledtab)[[1]])
	{
		for(b in dimnames(scaledtab)[[2]])
		{
			a = batches==b & conditions==c
			winnersamples = c( winnersamples, sample( samplenames[a], scaledtab[c,b]))
		}
	}
	return(winnersamples)	
}


#samplenames = sampleannotation$code
#batch = sampleannotation$chip
#covariate = sampleannotation$covariate1
#shufflecovariates = c("N", "DP")
# shuffle samples whithin a batch for selected covariates
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




xadhocpvalueplot = function(realcombatp, reallimmap, randomp, main="P-values", xrange=1:25)
{
  thiscolors = c("red", "blue", "black")
  
  a = hist(realcombatp, breaks=100, plot=F)$counts[xrange]
  b = hist(reallimmap, breaks=100, plot=F)$counts[xrange]
  c = hist(randomp, breaks=100, plot=F)$counts[xrange]
  ylim=c(0, max(c(a,b,c)))
  # reproduced ComBat + limma
  plot((xrange)/100, a , ylim=ylim,
       main=main, xlab="p-value", ylab="frequency",
       type="l", lwd=2, col=thiscolors[1])
  
  # batch handled in limma
  lines((xrange)/100, b,
        col=thiscolors[2], lwd=2)
  
  # random ComBat + limma
  lines((xrange)/100, c,
        col=thiscolors[3], lwd=2)
  
  legend("topright", text.col=thiscolors,
         legend=c("Real data, ComBat adjusted",
                  "Real data, batch handled by Limma",
                  "Random data, ComBat adjusted")) 
}


