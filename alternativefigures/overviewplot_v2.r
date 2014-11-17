
### Figure of box-plots and confidence intervals

library(limma)
library(lsmeans)

source("commonscripts/boxplot_function_v2.r")
source("commonscripts/helperfunctions.r")

# 3 batches.
# 4 groups
#sampleannotation = createsampleannotation(  list(c(100,0,0,20), c(0,100,0,20), c(0,0,100,20))) 
sampleannotation = createsampleannotation(  list(c(50,0,0,20), c(0,50,0,20), c(0,0,50,20))) 
#sampleannotation = createsampleannotation(  list(c(25,0,0,5), c(0,25,0,5), c(0,0,25,5))) 

table(sampleannotation[,2:3])
group=factor(sampleannotation$group)
batch=factor(sampleannotation$batch)

ngenes=1000
#matrix_random = createrandomdata(ngenes, sampleannotation, mean=0, sd=1)
rseed=1111
set.seed(rseed)


matrix_random = matrix(rnorm(ngenes * nrow(sampleannotation), mean=0, sd=1), nrow=ngenes, ncol=nrow(sampleannotation))

matrix_condition = matrix_random
matrix_condition[, group %in% c(3,4)] = matrix_condition[, group %in% c(3,4)] +3
matrix_conditionbatch=matrix_condition
matrix_conditionbatch[, batch ==2] = matrix_condition[, batch ==2] -1
matrix_conditionbatch[, batch ==3] = matrix_condition[, batch ==3] +1
matrix_meancenter = matrix_conditionbatch
matrix_meancenter[,batch ==1] =  matrix_meancenter[,batch ==1] - rowMeans( matrix_meancenter[,batch ==1])
matrix_meancenter[,batch ==2] =  matrix_meancenter[,batch ==2] - rowMeans( matrix_meancenter[,batch ==2])
matrix_meancenter[,batch ==3] =  matrix_meancenter[,batch ==3] - rowMeans( matrix_meancenter[,batch ==3])
mod = model.matrix(~group)
matrix_limma = removeBatchEffect(matrix_conditionbatch, batch=batch, design=mod)


index=1
#fit_lm=lm(matrix_condition[index,] ~ group+batch)
#means_lm=lsmeans(fit_lm,~group)
#summary(means_lm)
#contrast(means_lm,"pairwise")

figfilename = file.path( getwd(), "plots", paste("boxplots_v2", sep=""))
#figfile = paste( figfilename, ".png", sep=""); png(file = figfile, width=1600, height=800)
figfile = paste( figfilename, ".pdf", sep=""); 
pdf(file =figfile, width=24, height=24)

op=par(mfrow=c(5, 1), 	xaxt="n")
alldata = c(matrix_condition[index,],matrix_conditionbatch[index,],matrix_meancenter[index,],matrix_limma[index,])
ylim = c(min(alldata), max(alldata))
adhocboxplot(matrix_condition[index,], group, batch, ylim=ylim, xlab="Without batch effects", figureletter="a")
adhocboxplot(matrix_conditionbatch[index,], group, batch, ylim=ylim, xlab="With batch effects added", figureletter="b")
adhocboxplot(matrix_meancenter[index,], group, batch, ylim=ylim, xlab="Zero-centered per batch", figureletter="c")
adhocboxplot(matrix_limma[index,], group, batch, ylim=ylim, xlab="ANOVA centered values", figureletter="d")
adhocboxplot(matrix_conditionbatch[index,], group, batch, ylim=ylim, xlab="With batch effects added, ANOVA estimates", figureletter="e", estimatemethod="lsmeans")
#adhoclegend(group)
par(op)
dev.off()
print( paste("Figure created: ",figfile ))

