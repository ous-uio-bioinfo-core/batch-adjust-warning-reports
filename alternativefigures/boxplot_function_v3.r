



# y= matrix_true[index,]
# group=sampleannotation$group
# main=paste("True values")
# estimatemethod="CI"
# 
# # y=matrix_true[1,]
# # group=sampleannotation$group
# # batch=sampleannotation$batch
# batch=NULL
# # estimatemethod="lsmeans"
# ylim=NULL

plot_one_gene = function(y, group, batch=NULL, ylim=NULL, main="", estimatemethod="none", lwd=1)
{  
	adhoc.cex=1
	adhoc.legendcex=adhoc.cex
	adhoc.palette =  c("black","black","black","black")
	#adhoc.palette =c("red", "blue", "brown", "cyan")
	adhoc.pch = c(1, 2, 3, 4)
	
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
  xlim = c(1, round(xboxplots + (boxseparation) * length(unique(group)))  )
  
  # measurements
  plot(y, ylim=ylim, xlim=xlim, col=adhoc.palette[as.factor(group)],  pch=adhoc.pch[as.factor(group)],
       main=main,  ylab="Expression", lwd=lwd, cex.main=adhoc.legendcex)
  
  # batch separator
  if(!is.null(batch))
  {
    batchnames=unique(batch)
    x = match(batchnames, batch)
    abline(v = x[-1]-0.5, lty=3)#
    batchlabelx=vector()
    for(i in 1:(length(x)-1))
    {
      batchlabelx = c(batchlabelx, (x[i+1]+x[i])/2)
    }
    batchlabelx = c(batchlabelx, (x[i+1]+length(y))/2)	  
    text(labels="Batch", y=ylim[1]+0.5, x=2, cex=adhoc.legendcex)
    text(labels=paste(batchnames, sep=""), y=ylim[1]+0.5, x=batchlabelx, cex=adhoc.legendcex)	  
  }
  
  # group means line
  groupnames=factor(unique(group))
  groupmeans = vector()
  for(thisgroup in groupnames)
  {
    groupmeans[thisgroup]=round(mean(y[group==thisgroup]),2)
  }  
  segments(  x0=(1:length(y))-0.5, y0=groupmeans[group], 
             x1=(1:length(y))+0.5, y1=groupmeans[group], 
             col=adhoc.palette[as.factor(group)], lwd=lwd, lty=2)
  #legend("topright", legend=paste(groupnames, groupmeans), text.col=adhoc.palette[as.factor(groupnames)])
  
  
  
  #text(x=xoffset+boxwidth, y=ylim[2]-0.5, labels=estimatemethod, cex=legendcex)
  xoffset= xboxplots
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
    xleft = xoffset    
    #print(xleft)
    xright = xleft + boxwidth    
    rect(xleft, ybottom, xright, ytop, border=adhoc.palette[g], lwd=lwd, lty="solid")
    lines(c(xoffset,xright), c(m,m), col=adhoc.palette[g], lwd=lwd)
    points(x=xleft+boxwidth/2, y=ylim[1], pch=adhoc.pch[g], cex=adhoc.legendcex, lwd=lwd)
    xoffset = xoffset+boxseparation
  }
  
  if(estimatemethod=="lsmeans")
  {
    pval= round(summary(contrast(means_lm,"pairwise"))$p.value[1],2)
    fc= round(summary(contrast(means_lm,"pairwise"))$estimate[1],2)    
  } else if(estimatemethod=="CI") {
    fc = round(groupmeans[1]-groupmeans[2],2)
    pval= round(t.test( y[group==groupnames[1]], y[group==groupnames[2]])$p.value,2)
  } 
  difftestlabel = paste("FC=", fc, ", p=", pval, sep="")
  text(labels=difftestlabel, y=ylim[2]-0.5, x=xboxplots+boxwidth, cex=adhoc.legendcex) 
  #legend("topright", legend=difftestlabel, cex=adhoc.legendcex) 
  
  
  
}
