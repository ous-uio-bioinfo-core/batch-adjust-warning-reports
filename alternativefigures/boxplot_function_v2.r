# Functions for producing a boxplot

#adhoc.palette =  c("blue", "red", "darkgreen", "darkorange")
adhoc.palette =  c("black","black","black","black")
adhoc.pch = c(1, 2, 3, 4)
adhoc.cex=3
adhoc.legendcex=adhoc.cex
adhoc.labcex=adhoc.cex

# y=matrix_conditionbatch[index,]
# grouplabels=group
# batchlabels=batch
# addlegend=FALSE
# ylim=ylim
# xlab="Without batch effects"
# figureletter="a"
# cex=adhoc.cex
# labelcex=adhoc.labcex
# legendcex=adhoc.legendcex
# estimatemethod="lsmeans"

# Boxplot component: values and intervals
adhocboxplot = function(y, grouplabels, batchlabels, addlegend=FALSE, ylim=NULL, xlab="x",figureletter=NA,cex=adhoc.cex,labelcex=adhoc.labcex, legendcex=adhoc.legendcex, estimatemethod="CI")
{
  colpalette = adhoc.palette
  pchselection = adhoc.pch

  lwd=2
  
  xlim = c(1, round(length(y)*1.25))
  plot(y, col=colpalette[grouplabels], main=xlab, cex.main=labelcex, pch=pchselection[grouplabels],
       xlim=xlim, lwd=lwd, cex=cex, ylim=ylim, ylab="", xlab=NA, cex.lab =labelcex)

  #add batch separators
  for(thisbatch in unique(batchlabels))
  {
    sepx = length(batchlabels) - match(thisbatch, batchlabels[length(batchlabels):1])+1+0.5
    abline(v=sepx, col="black", lwd=lwd)
    text(x=sepx-length(batchlabels)/7, y=ylim[2]-0.5, labels=paste("Batch ", thisbatch, sep=""), cex=legendcex)    
  }
  
  #add boxplot
  groupnames = sort(unique(grouplabels))
  xoffset=round(length(y)* 1.1)  
  boxwidth=round(length(y)/20)
  boxseparation = round(boxwidth/1.5)
  text(x=xoffset+boxwidth, y=ylim[2]-0.5, labels=estimatemethod, cex=legendcex)
  for(g in groupnames) {
    g=factor(g,levels=levels(groupnames))
    
    if(estimatemethod=="lsmeans")
    {
      fit_lm=lm(y ~ grouplabels+batchlabels)
      means_lm=lsmeans(fit_lm,~grouplabels)
      anovaest = summary(means_lm)
      g=factor(g,levels=levels(groupnames))
      m =anovaest[g , "lsmean"]      
      ybottom = anovaest[g , "lower.CL"]      
      ytop = anovaest[g , "upper.CL"]      
    }
    else
    {
      m = mean(y[grouplabels==g])
      ci = t.test(y[grouplabels==g])[["conf.int"]]
      ybottom = ci[1]
      ytop = ci[2]
    }
    
    xleft = xoffset    
    xright = xleft + boxwidth    
    rect(xleft, ybottom, xright, ytop, border=colpalette[g], lwd=lwd, lty="solid")
    lines(c(xoffset,xright), c(m,m), col=colpalette[g], lwd=lwd)
    points(x=xleft+boxwidth/2, y=ylim[1]+0.5, pch=pchselection[g], cex=legendcex, lwd=lwd)
    xoffset = xoffset+boxseparation
  }
  
  
	if(!is.na(figureletter)) {
		abcde=paste("(", figureletter, ")", sep="")
		text(x=1, y=ylim[2]-0.5, labels=abcde, cex=legendcex)
	}
}



adhoclegend=function(groups,legendcex=adhoc.legendcex) {
	groupnames=sort(unique(groups))
	legend("topright", legend=paste("Group ", groupnames, sep=""), 
				 text.col=adhoc.palette[groupnames], bty="n", cex=legendcex)

}
