library(Biobase)
library(sva)
library(limma)
library(lsmeans)
library(pvca)
library(dendextend)

makeoptionslink = function(input)
{
	inputlist = reactiveValuesToList(input)
	link = "?"
	for(n in names(inputlist)[order(names(inputlist))])
		link = paste(link, n, "=", inputlist[[n]][1], "&", sep="")
	
	ret = paste("[ <a href='",link,"'>Link to current parameter settings</a> ]", sep="")
	return(ret)
}

	
batchadjustplot = function(input)
{ 
  print("i batchadjustplot")
  set.seed(input$rngseed)
  #ngenes = 10000
  ngenes = input$ngenes
  #index=input$indexgene
  index=1
  
  inputlist = reactiveValuesToList(input)

  combocounts=inputlist[grepl("count", names(inputlist))]
  names(combocounts)= gsub("count", "", names(combocounts))
  sa = createsampleannotation(  combocounts, as.factors=TRUE)
  rownames(sa) = sa$id
  matrix_true = matrix(rnorm(ngenes * nrow(sa), mean=10, sd=2), nrow=ngenes, ncol=nrow(sa))
  colnames(matrix_true) = rownames(sa)
  
  # add a group effect
  groupeffects = inputlist[grepl("^groupeffect[ABC]$", names(inputlist))]
  names(groupeffects) = gsub("groupeffect", "",names(groupeffects) )
  DEGindexes = 1:round(nrow(matrix_true) * as.numeric(input$groupeffectfraction))
  for(group in names(groupeffects))
  {   
    thiseffects = groupeffects[[group]]
    if(input$groupeffecttype=="sd")
    {
      thiseffects = rnorm(length(DEGindexes), mean=0, sd=abs(groupeffects[[group]]))
    }    
    matrix_true[DEGindexes,sa$group==group] = matrix_true[DEGindexes,sa$group==group]+thiseffects
  }
  if(input$zerocentre)
  	matrix_true = matrix_true - rowMeans(matrix_true)

  
  # calculating batchbox limits.
  bbextra=1.2
  batchboxheight=(max(matrix_true[index,]) - min(matrix_true[index,]) ) * bbextra
  batchboxlow=min(matrix_true[index,]) - batchboxheight * (bbextra-1)/2
  batchboxtop=max(matrix_true[index,]) + batchboxheight * (bbextra-1)/2
  batchmeans = sapply( unique(sa$batch), 
                       FUN=function(x){  mean( matrix_true[index,sa$batch==x] ) } )
  batchboxlowmeanoffsets=batchmeans - batchboxlow
  
  cat(input$adjustmethod)
  
  # add batch effects
  batcheffects = inputlist[grepl("batcheffectvalue", names(inputlist))]
  batcheffectsadd = matrix_true
  batcheffectsadd[]=0
  for(be in names(batcheffects))
  {
    thiseffects = batcheffects[[be]]
    if(input$batcheffecttype=="sd")
    {
    	thiseffects = rnorm(nrow(batcheffectsadd), mean=0, sd=abs(batcheffects[[be]]))
    }
    	
    batch = gsub("batcheffectvalue", "", be)
    batcheffectsadd[,sa$batch==batch] = thiseffects
    #batcheffectsadd[,sa$batch==batch] =rnorm(length(batcheffectsadd[,sa$batch==batch]), mean=batcheffects[[be]], sd=3^as.numeric(batch))# ad hoc test of different variation in the batches.
  }
  matrix_batcheffect = matrix_true + batcheffectsadd
  if(input$zerocentre)
  	matrix_batcheffect = matrix_batcheffect - rowMeans(matrix_batcheffect)
  

  matrix_batchadjusted = matrix_batcheffect
  if(input$adjustmethod=="Mean-centring")
  {
    # batch adjust with mean center.    
    genemeans = rowMeans(matrix_batcheffect)
    for(batch in unique(sa$batch))
    {
      matrix_batchadjusted[,sa$batch ==batch] =  genemeans +
        matrix_batchadjusted[,sa$batch ==batch] - rowMeans( matrix_batchadjusted[,sa$batch ==batch])
    }
  }
  if(input$adjustmethod=="removeBatchEffect")
  {
    # batch adjust with removeBatchEffect
  	mod = model.matrix(~as.factor(sa$group))
    matrix_batchadjusted = removeBatchEffect(matrix_batcheffect, batch=as.factor(sa$batch), design=mod)
    
  }
  if(input$adjustmethod %in% c("ComBat", "ComBat no covariates"))
  {
    # batch adjust with ComBat
  	combatmod = model.matrix(~as.factor(sa$group))
  	if(input$adjustmethod =="ComBat no covariates")
  		combatmod = NULL
    matrix_batchadjusted = ComBat(dat=matrix_batcheffect, batch=as.factor(sa$batch), 
                                  mod=combatmod, numCovs=NULL, 
                                  par.prior=TRUE, prior.plots=FALSE)    
  }
  if(input$zerocentre)
  	matrix_batchadjusted = matrix_batchadjusted - rowMeans(matrix_batchadjusted)
 
  pair=c()
  tpindexmax=0
  if(input$plotpvaluehist %in% c("AB", "AC", "BC"))
  {
  	pair = strsplit(input$plotpvaluehist, "")[[1]][1:2]
		if(groupeffects[[pair[1]]]!=groupeffects[[pair[2]]])
			tpindexmax = max(DEGindexes)
  }
 
  
  #op=par(mfrow=c(3, 2))
  nrows = as.numeric(input$plottrue)+as.numeric(input$plotbatchaffected)+as.numeric(input$plotbatchadjusted)+as.numeric(input$plotlsmeans)
  ncols = 1+as.numeric(length(pair)>0)
  nplots = nrows * ncols
  
  if(length(pair)>0)
    widths=c(1,0.35)
  else
    widths=c(1)
  
  
  #layout(  matrix( 1:(cols*rows), rows, cols, byrow = TRUE) , widths=widths )
  if(input$blackbg)
    par(bg= "black", fg="white", col.axis="white", col.lab="white", col.main="white", col.sub="white")
  else
    par(bg= "white", fg="black", col.axis="black", col.lab="black", col.main="black", col.sub="black")
  
   save(matrix_true, matrix_batcheffect, matrix_batchadjusted, sa, batchboxheight, batchboxlowmeanoffsets, input, index, file="not_in_github/image.rdata")



batchmean_real=unlist(lapply(unique(sa$batch), FUN=function(x) mean(matrix_true[index,sa$batch==x])))
batchmean_batch =unlist(lapply(unique(sa$batch), FUN=function(x) mean(matrix_batcheffect[index,sa$batch==x])))
batchmean_adjusted =unlist(lapply(unique(sa$batch), FUN=function(x) mean(matrix_batchadjusted[index,sa$batch==x])))  

  layout(  matrix( 1:8, 4, 2, byrow = TRUE) , widths=c(1,0.35) )
  shouldaddlegend=TRUE
  if(input$plottrue)
  {
    plot_one_gene(matrix_true[index,], group=sa$group,  batch=sa$batch,, 
                main=paste("\"True\" values, one gene", sep=""), estimatemethod="CI", ylim=NULL,
                bbh=batchboxheight, bblo = batchboxlowmeanoffsets, addlegend=shouldaddlegend,
    							plotCI=input$plotCI, addbatchbox=input$plotbatchbox)
    shouldaddlegend=FALSE
    if(length(pair)>0)
      plot_pvals(matrix_true, group=sa$group, pair=pair, pvaluemax=input$pvalueplotxlim, tpindexmax=tpindexmax)
    else if(input$plotpvaluehist=="pvca")
      plot_pvca(matrix_true, sa[,-1])
    else if(input$plotpvaluehist=="pca")
    	plot_pca(matrix_true, sa[,-1])
    else if(input$plotpvaluehist=="hclust")
      plot_hclust(matrix_true, sa[,-1]) 
    else
      plot_blank()
  }
  if(input$plotbatchaffected)
  {
  	bbshift=NULL
  	if(input$printbatchshift)
  		bbshift=batchmean_batch-batchmean_real
    plot_one_gene(matrix_batcheffect[index,], group=sa$group, batch=sa$batch,
                main=paste("Observed values (batch affected), one gene", sep=""), estimatemethod="CI", 
                bbh=batchboxheight, bblo = batchboxlowmeanoffsets, bbshift=bbshift,addlegend=shouldaddlegend,
    							plotCI=input$plotCI, addbatchbox=input$plotbatchbox)
    shouldaddlegend=FALSE
    if(length(pair)>0)
      plot_pvals(matrix_batcheffect, group=sa$group, pair=pair, pvaluemax=input$pvalueplotxlim, tpindexmax=tpindexmax)
  	else if(input$plotpvaluehist=="pvca")
  	  plot_pvca(matrix_batcheffect, sa[,-1])
  	else if(input$plotpvaluehist=="pca")
  		plot_pca(matrix_batcheffect, sa[,-1])
  	else if(input$plotpvaluehist=="hclust")
  	  plot_hclust(matrix_batcheffect, sa[,-1]) 
  	else
  	  plot_blank()
  }
  
  if(input$plotbatchadjusted)
  {
  	bbshift=NULL
  	if(input$printbatchshift)
  		bbshift=batchmean_adjusted-batchmean_batch
    plot_one_gene(matrix_batchadjusted[index,], group=sa$group, batch=sa$batch,
                main=paste(input$adjustmethod, " adjusted values, one gene", sep=""), estimatemethod="CI",
                bbh=batchboxheight, bblo = batchboxlowmeanoffsets, bbshift=bbshift,addlegend=shouldaddlegend,
    							plotCI=input$plotCI, addbatchbox=input$plotbatchbox)
    shouldaddlegend=FALSE
    if(length(pair)>0)
      plot_pvals(matrix_batchadjusted, group=sa$group, pair=pair, pvaluemax=input$pvalueplotxlim, tpindexmax=tpindexmax)
  	else if(input$plotpvaluehist=="pvca")
  	  plot_pvca(matrix_batchadjusted, sa[,-1])
  	else if(input$plotpvaluehist=="pca")
  		plot_pca(matrix_batchadjusted, sa[,-1])
  	else if(input$plotpvaluehist=="hclust")
  	  plot_hclust(matrix_batchadjusted, sa[,-1])  	
  	else
  	  plot_blank()
  }

  if(input$plotlsmeans)
  {
    # lsmeans 2way anova
    estimatesboxesonly(matrix_batcheffect[index,], group=sa$group, batch=sa$batch, 
                       main="Batch as factor in two way ANOVA. LSmeans estimates 95% CI")
    shouldaddlegend=FALSE

    if(length(pair)>0)
      plot_pvals(matrix_batcheffect, group=sa$group, batch=sa$batch, 
                pair=pair, pvaluemax=input$pvalueplotxlim, tpindexmax=tpindexmax)
    else
      plot_blank()
  }
  
  
  
  
  nblankplots = (4 - nrows) * 2
  for (i in seq_along(  numeric(nblankplots) ))
    plot_blank()

}



##### plot visual settings
adhoc.cex=2
adhoc.legendcex=adhoc.cex
adhoc.pch = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3)

adhoc.usecolor = function(usecolor=TRUE) {
	if (usecolor) {
		adhoc.palette <<- c("darkred", "darkgreen", "darkblue")
		adhoc.palette.line <<- c("red", "seagreen", "blue"); # adjustcolor(adhoc.palette[g], linealpha)
		adhoc.palette.fill <<- c("lightpink", "seagreen1", "lightblue"); # adjustcolor(adhoc.palette[g],boxalpha))
		adhoc.batch.colour<<- "darkgray"; #adjustcolor("black", 0.3)
	} else {
		adhoc.palette <<-  c("black","black","black","black")
		adhoc.palette.line <<- c("darkgray", "darkgray", "darkgray"); # adjustcolor(adhoc.palette[g], linealpha)
		adhoc.palette.fill <<- c("lightgray", "lightgray", "lightgray"); # adjustcolor(adhoc.palette[g],boxalpha))
		adhoc.batch.colour <<- "darkgray"; #adjustcolor("black", 0.3)
	}
}
adhoc.usecolor();



plot_blank = function()
{
  plot(0,type='n',axes=FALSE,ann=FALSE)
}

plot_pvals = function(dm, group, pair, batch=NULL, pvaluemax=1, fdr=0.05, tpindexmax=0)
{
  pvaluemax = as.numeric(pvaluemax)
  group = factor(group)
  design = model.matrix(~0 + group)
  if(!is.null(batch))
    design = model.matrix(~0+group+batch)
  fit = lmFit(dm, design)
  contrast = paste( "group", pair[1], "-", "group", pair[2], sep="")
  cont.matrix = makeContrasts ( contrasts=contrast, levels=design)  
  fit2 = contrasts.fit(fit, cont.matrix)
  limma_p = eBayes(fit2)$p.value[,1]
  limma_fdr = p.adjust(limma_p, method="BH")
  deg = sum(limma_fdr<fdr)
  tp = 0
  if(tpindexmax>0)
 		tp = sum(limma_fdr[1:tpindexmax]<fdr)
  
  breaks = 20 * 1/pvaluemax
  hist(limma_p, breaks=breaks, xlim=c(0,pvaluemax), 
       main=paste(pair[1], "vs." , pair[2], ", ", nrow(dm), " simulations", sep=" "), xlab="P-value")
  legendstrings = c( paste("FDR<", fdr, ": ", deg, sep=""),
  									 paste("True Positive: ", tp, sep="")
  									 )
  legend("topright", legendstrings, bty="n")
  
}


plot_pvca = function(dm, sa)
{
 
  pct_threshold=0.6
  eset <- ExpressionSet(assayData=dm, phenoData=new("AnnotatedDataFrame", data=sa))
  pvcaObj = pvcaBatchAssess (eset, c("batch", "group"), pct_threshold) 
 
  bp <- barplot(pvcaObj$dat,  xlab = "", 
                ylab = "Weighted average proportion variance", 
                ylim= c(0,1.1),col = c("grey"), las=2, 
                main="PVCA estimation")
  axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 1.5, las=2)
  values = pvcaObj$dat
  new_values = round(values , 3)
  text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 1) 
}

plot_pca = function(dm, sa)
{
	thisprcomp=prcomp( t(dm) )
	plot(thisprcomp$x,  col=adhoc.palette[sa$group], cex=2,  main="PCA-plot", pch=as.character(sa$batch) , cex.main=1, axes=FALSE)
}

plot_hclust = function(dm, sa)
{	
	dend <- as.dendrogram(hclust(dist(t(dm))))
	groupcols = adhoc.palette[ sa[labels(dend), "group"]]
	dend <- color_branches(dend, col =groupcols )
	if(ncol(dm)>25)
	{
		labels(dend) = paste( " ", unlist(lapply( sa[labels(dend), "batch"], FUN=function(x)paste(rep("-", x), collapse="") )), sep="")
	}	else {
		labels(dend) = as.character(sa[labels(dend), "batch"])
		dend <- set(dend, "leaves_pch", 15)
		dend <- set(dend, "leaves_cex", 3)
	}
	labels_colors(dend) = "black"
	dend <- set(dend, "labels_cex", 1.5)
	dend <- set(dend, "branches_lwd", 3)
	dend <- set(dend, "leaves_col", groupcols )
	plot(dend, main = "", axes=F, horiz =  TRUE)
}

plot_one_gene = function(y, group, batch=NULL, ylim=NULL, main="", 
                         estimatemethod="none", lwd=1, boxlabel='CI',leftmargin=1,
                         bbh=NA, bblo=NA, bbshift=NULL, addlegend=FALSE, plotCI="both", addbatchbox=FALSE)
{  
  
  
  # Boxplots for CI etc.  
  xboxplots=round(length(y)* 1.1)    
  boxwidth=round(length(y)/20)
  boxseparation = (boxwidth * 1.5)
  
  
  a = order(group)
  if(!is.null(batch))
  {
    a = order(batch, group)
    batch=batch[a]
  }
  group=group[a]
  
  y=y[a]  
  
  if(is.null(ylim))
  {
    ymax = round(max(y)+2, 1)
    ymin = round(min(y)-2, 1)
    ylim=c(ymin,ymax)
  }
  #xlim = c(1, round(xboxplots + (boxseparation) * length(unique(group)))  )
  
  rightboxwidth = (length(y)/6)/length(unique(group))
  rightboxsep=rightboxwidth/5
  rightboxsoffset = 1.5 * rightboxwidth
  rgihtboxtotal = rightboxsoffset+ (rightboxwidth + rightboxsep)*length(unique(group))
  xlim = c(1-leftmargin, length(y)*1.03)
  if(plotCI %in% c("both", "right"))
    xlim[2] = xlim[2]+rgihtboxtotal
  
    
  
  # measurements
  plot(y, ylim=ylim, xlim=xlim, col=adhoc.palette[as.factor(group)],  pch=adhoc.pch[as.factor(group)], main=main,  ylab=NA, xlab=NA, lwd=lwd, cex.main=adhoc.legendcex, xaxt="n", yaxt="s", cex=adhoc.cex, type="n")
  
  
  
  
  groupnames=factor(unique(group))
  # color code legend
  if(addlegend)
  {
    legend("topright", legend=groupnames, text.col=adhoc.palette[groupnames], cex=adhoc.cex, bty="n")
  }
  
  currentrightboxx=length(y) + rightboxsoffset 
  for(g in groupnames)
  {
    g=factor(g,levels=levels(groupnames))
    
    if(estimatemethod=="lsmeans")
    {
      fit_lm=lm(y ~ group+batch)
      means_lm=lsmeans(fit_lm,~group)
      anovaest = summary(means_lm)      
      m =anovaest[g , "lsmean"]      
      ybottom = anovaest[g , "lower.CL"]      
      ytop = anovaest[g , "upper.CL"]      
    } else if(estimatemethod=="CI") {
      m = mean(y[group==g])
      ci = t.test(y[group==g])[["conf.int"]]
      ybottom = ci[1]
      ytop = ci[2]
    }
    
    if(plotCI %in% c("both", "over"))
    {
      xextra = 0.3
      xleft = match(g, group) - xextra
      xright = length(group) - match(g, group[length(group):1]) + 1 + xextra
    
      linealpha=0.2
      boxalpha=0.1
      rect(xleft, ybottom, xright, ytop, border=adhoc.palette.line[g], lwd=lwd, 
           lty="solid", density=-1, col=adhoc.palette.fill[g])
      lines(c(xleft,xright), c(m,m), col=adhoc.palette.line[g], lwd=lwd)
    }
    if(plotCI %in% c("both", "right"))
    {
      xleft = currentrightboxx
      xright = xleft + rightboxwidth
      rect(xleft, ybottom, xright, ytop, border=adhoc.palette[g], lwd=lwd*2, 
           lty="solid", density=-1, col=adhoc.palette.fill[g])
      
      lines(c(xleft,xleft,xright,xright,xleft), c(ytop,ybottom,ybottom,ytop,ytop), col=adhoc.palette[g], lwd=lwd*2)
      lines(c(xleft,xright), c(m,m), col=adhoc.palette[g], lwd=lwd*2)
      currentrightboxx = currentrightboxx + rightboxwidth + rightboxsep
    }
    
    
    # batch separator
    if(!is.null(batch) & addbatchbox)    
    {
    	for(b in unique(batch))
    	{
    		batchnames=unique(batch)
    		#x = match(batchnames, batch)
    		#abline(v = x[-1]-0.5, lty=3)#
    		
    		xleft = match(b, batch)-0.4
    		xright = length(batch) - match(b, batch[length(batch):1]) + 1 +0.4
    		batchmean = mean(y[batch==b])
    		ybottom = batchmean - bblo[as.numeric(b)]
    		ytop = ybottom + bbh    
    		rect(xleft, ybottom, xright, ytop, lty=3, border=adhoc.batch.colour)
  		
				# Decide where to put batch label, up of down.  Try to print inside the ylim.
    		if( ybottom > ylim[1] ) # bottom
     		{     			
     			laby=ybottom
     		}else{ # top
     			laby=ytop
     		}

    		shiftstring = ""
    		if(!is.null(bbshift))
    		{
    			if(bbshift[as.numeric(b)] > 0)
    				shiftstring=paste( "+" , round(bbshift[as.numeric(b)], 2), sep="")
    			else
    				shiftstring=paste(  round(bbshift[as.numeric(b)],2), sep="")
    		}      
    		text(labels=paste("  Batch", b,"  ",shiftstring , sep=""), y=laby, x=xleft, pos = 4 ,
    				 cex=adhoc.legendcex, col="black")
    	}
    }
    
    points(y, col=adhoc.palette[as.factor(group)],  pch=adhoc.pch[as.factor(group)], lwd=2, cex=3)

    
  } 
  
}



estimatesboxesonly = function(y, group, batch, ylim=NULL, main="", lwd=1)
{
  
  groupnames=factor(unique(group))
  
  if(is.null(ylim))
  {
    ymax = round(max(y)+2, 1)
    ymin = round(min(y)-2, 1)
    ylim=c(ymin,ymax)
  }
  
  plot( as.numeric(unique(group)),  xlim=c(1, length(unique(group))+1) , 
        ylim=ylim, type="n", xaxt="s", yaxt="s", main=main, cex.main=adhoc.legendcex, xlab=NA, ylab=NA)
  for(g in groupnames)
  {
    g=factor(g,levels=levels(groupnames))
    
    
    fit_lm=lm(y ~ group+batch)
    means_lm=lsmeans(fit_lm,~group)
    anovaest = summary(means_lm)      
    m =anovaest[g , "lsmean"]      
    ybottom = anovaest[g , "lower.CL"]      
    ytop = anovaest[g , "upper.CL"]      
    
    
    xleft = as.numeric(g)
    xright = as.numeric(g) + 0.9
    
    linealpha=0.2
    boxalpha=0.1
    rect(xleft, ybottom, xright, ytop, border=adhoc.palette.line[g], lwd=lwd, 
         lty="solid", density=-1, col=adhoc.palette.fill[g])
    lines(c(xleft,xright), c(m,m), col=adhoc.palette.line[g], lwd=lwd)
    #points(x=(xleft+xright)/2, y=ytop+( (ylim[2]-ylim[1])/30), pch=adhoc.pch[g], cex=adhoc.legendcex, lwd=lwd, col=adhoc.palette[g])
    
  }
}


#groupbatchcount = list(A2=0, A3=0, B1=2, B2=5, B3=2, C1=0, C2=0, C3=4, A1=3)
#groupbatchcount = list(A1=5, A2=3, B1=2, B2=4)

# Creates a sampleannotation data.frame based on a list of counts where the name is the batchgroup combo.
# for instance  A1=5, A2=3, B1=3, B2=4
# Used when creating artificial data.
createsampleannotation = function( groupbatchcount, as.factors=FALSE)
{
  n<-sum(unlist(groupbatchcount))
  sampleannotation = data.frame( id=paste("sample", 1:n, sep=""),    group = character(n), batch=numeric(n), stringsAsFactors=FALSE)
  
  count<-1
  for(combo in names(groupbatchcount))
  {
    if(groupbatchcount[[combo]][[1]] > 0 )
    {
      start<-count
      end<-start+groupbatchcount[[combo]][[1]]-1    
      #sampleannotation[start:end, "group"] = strsplit(combo, "")[[1]][1]
      #sampleannotation[start:end, "batch"] = strsplit(combo, "")[[1]][2]
      sampleannotation[start:end, "group"] = substring(combo, 1,1)
      sampleannotation[start:end, "batch"] = substring(combo, 2)
      count<-end+1
    }  
  }
  
  if(as.factors)
  {
    sampleannotation$batch=as.factor(sampleannotation$batch)
    sampleannotation$group=as.factor(sampleannotation$group)
  }
  sampleannotation = sampleannotation[order(sampleannotation$batch, sampleannotation$group),]
  return(sampleannotation)
}
