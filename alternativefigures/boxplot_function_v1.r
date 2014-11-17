# Functions for producing a boxplot

adhoc.palette =  c("blue", "red", "darkgreen", "darkorange")
adhoc.pch = c(1, 2, 3)
adhoc.cex=3
adhoc.legendcex=adhoc.cex
adhoc.labcex=adhoc.cex

# Boxplot component: values and intervals
adhocboxplot = function(y, grouplabels, batchlabels, addlegend=FALSE, ylim=NULL, xlab="x",figureletter=NA,cex=adhoc.cex,labelcex=adhoc.labcex, legendcex=adhoc.legendcex)
{
  colpalette = adhoc.palette
  pchselection = adhoc.pch
  x = as.numeric(batchlabels)
  lwd=2
  
  # add extra x for better visual separation
  for(b in unique(batchlabels)) {
    separation = 0.15
    xoffset=0
    for(g in unique(grouplabels[batchlabels==b])) {
      x[batchlabels==b & grouplabels==g] = x[batchlabels==b & grouplabels==g] + xoffset
      xoffset = xoffset+separation
    }
  }
  
  xlim = c(0, length(unique(batchlabels))+4)
  plot(x, y, col=colpalette[grouplabels], pch=pchselection[batchlabels], 
       xlim=xlim, lwd=lwd, cex=cex, ylim=ylim, xlab=xlab, ylab="", cex.lab =labelcex)

  #add boxplot
  groupnames = sort(unique(grouplabels))
  xoffset=length(unique(batchlabels))+2
  boxwidth=0.35
  boxseparation = 0.25
  for(g in groupnames) {
    g=factor(g,levels=levels(groupnames))
    m = mean(y[grouplabels==g])
    ci = t.test(y[grouplabels==g])[["conf.int"]]
    xleft = xoffset
    ybottom = ci[1]
    xright = xleft + boxwidth
    ytop = ci[2]
    rect(xleft, ybottom, xright, ytop, border=colpalette[g], lwd=lwd)
    lines(c(xoffset,xright), c(m,m), col=colpalette[g], lwd=lwd)
    xoffset = xoffset+boxseparation
  }
	if(!is.na(figureletter)) {
		figureletter=paste("(", figureletter, ")", sep="")
		text(x=1, y=ylim[2], labels=figureletter, cex=legendcex)
	}
}


# Boxplot component: only intervals
adhocboxplot2 = function(lsmeans, grouplabels, batchlabels, addlegend=FALSE, ylim=NULL, xlab="x", figureletter=NA, cex=adhoc.cex,labelcex=adhoc.labcex,legendcex=adhoc.legendcex)
{
	colpalette = adhoc.palette
	pchselection = adhoc.pch
	x = as.numeric(batchlabels)
	lwd=2
	
	xlim = c(0, length(unique(batchlabels))+4)
	plot(1, 1, col=colpalette[grouplabels], pch=pchselection[batchlabels], 
			 xlim=xlim, lwd=lwd, cex=cex, ylim=ylim, xlab=xlab, ylab="", cex.lab=adhoc.labcex, type="n")
	
	anovaest = summary(lsmeans)
	
	groupnames = sort(unique(grouplabels))
	xoffset=2
	boxwidth=0.35
	boxseparation = 0.25
	for(g in groupnames) { 
		g=factor(g,levels=levels(groupnames))
		m =anovaest[g , "lsmean"]
		xleft = xoffset
		ybottom = anovaest[g , "lower.CL"]
		xright = xleft + boxwidth
		ytop = anovaest[g , "upper.CL"]
		rect(xleft, ybottom, xright, ytop, border=colpalette[g], lwd=lwd)
		lines(c(xoffset,xright), c(m,m), col=colpalette[g], lwd=lwd)
		xoffset = xoffset+boxseparation
	}
	if(!is.na(figureletter)) {
		figureletter=paste("(", figureletter, ")", sep="")
		text(x=1, y=ylim[2], labels=figureletter, cex=legendcex)
	}
}

adhoclegend=function(groups,batches,legendcex=adhoc.legendcex) {
	groupnames=sort(unique(groups))
	batchnames=sort(unique(batches))
	legend("topright", legend=paste("Group ", groupnames, sep=""), 
				 text.col=adhoc.palette[groupnames], bty="n", cex=legendcex)
	legend("bottomright", legend=paste("Batch ", batchnames, sep=""),
				 pch=adhoc.pch[batchnames], bty="n", cex=legendcex)
}
