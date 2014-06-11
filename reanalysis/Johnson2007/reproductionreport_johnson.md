

Reproduction and alternative analysis of "Data Set 2" from "Adjusting batch effects in microarray data using Empirical Bayes methods."
========================================================



2014-06-11 14:13:43

### Overview
This report aims to show to what extent the use of ComBat led to false results in the second analysis example given in [Johnson et al.](http://biostatistics.oxfordjournals.org/content/8/1/118.abstract) The example named "Data Set 2" and the analysis is described in the [supplementary material](http://biostatistics.oxfordjournals.org/content/suppl/2006/04/21/kxj037.DC1/kxj037supp.pdf) for Johnson et al.
Description on how to obtain the data were found [here](https://groups.google.com/d/msg/combat-user-forum/S9vBcXw8RGk/XNWnNl0PE94J) which had links to [dataExample2.txt](http://www.bu.edu/jlab/wp-assets/ComBat/data/dataExample2.txt)
and [sampleInfoExample2.txt](http://www.bu.edu/jlab/wp-assets/ComBat/data/sampleInfoExample2.txt)

This document has four main parts 
- Reproduce some of the results to show that we are working on the same data and analysis workflow
- Remove the use of ComBat and perform the same analysis with an alternative established tool
- Estimate the error introduced when ComBat is used and the consequences for the conclusion of the study
- Perform a few more sanity checks to substantiate that the difference in results for the two above analyses is mainly false and introduced by ComBat  

        
### Read data and sample annotation
The data files consist of 35 samples of which 5 are annotated as "WT" and were not
referred to in Johnson et al. These are taken out in the beginning and are not used in this report.


```r
includelibs = c("pheatmap", "sva", "qvalue", "limma", "NMF", "statmod")
# NB! the qvalue library has a dependency (maybe tcl-related) that starts a
# x11 icon.  Closing that x11 icon will abort the R session. Otherwise it
# seem to work.
tmp = lapply(includelibs, require, character.only = T)
print(tmp)
if (any(!unlist(tmp))) {
    stop(paste("Not able to find all packages. Please install ", paste(includelibs[!unlist(tmp)], 
        collapse = ", ")))
}
rm(tmp)
debug = FALSE
```


```r
datamatrix = as.matrix(read.table("data/dataExample2.txt", sep="\t", header=TRUE))
if(debug)datamatrix = datamatrix[1:3000,]
sampleannotation = read.table("data/sampleInfoExample2.txt", 
                              sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(sampleannotation)=sampleannotation$ArrayName
sampleannotation$Batch=factor(as.character(sampleannotation$Batch)) 
sampleannotation$Cell=factor(as.character(sampleannotation$Cell)) 
# must be discrete for the pheatmap
#table(sampleannotation$ArrayName==dimnames(datamatrix)[[2]])#ordercheck
datamatrix=datamatrix[,sampleannotation$Type!="WT"]
sampleannotation=sampleannotation[sampleannotation$Type!="WT",]
sampleannotation$Type=factor(sampleannotation$Type)
#dev/debug
useparprior=TRUE
print(dim(datamatrix))
```

```
## [1] 54675    30
```

Inspection of the covariate/batch balance;

```r
print(table(sampleannotation[, c("Batch", "Type")]))
```

```
##      Type
## Batch C R
##     1 2 6
##     2 4 3
##     3 6 9
```

Three batches and two covariates. There is a "Cell" column in the sample annotation, but it does not seem to have been used.
<br/>  
<br/>  
### Reproduce the original results
Following the description in section A.1 and A.2 and figure texts in the supplementary materials for Johnson et al., we try to reproduce some of their results. First, the heatmap in Figure A.1.
> A heatmap clustering of data set 2. 698 genes with large variation across all the samples are clustered.<cite> Johnson et al.


```r
variationmeasure = apply(datamatrix, 1, FUN= function(x){var(x)})
clustermatrix = datamatrix[order(variationmeasure, decreasing=TRUE),][1:2698,]
aheatmap(clustermatrix, scale = "row", Rowv=T, Colv=T, 
         color = colorRampPalette(c("red", "black", "green"))(n = 299),
         annCol = sampleannotation[, c("Batch", "Type")], 
         cellheight=NA, cellwidth=13,fontsize=12,   labRow=NA,
         main=paste("Reproduced Fig A1.",sep=""),border_color=NA,treeheight=c(0,100) )

rm(variationmeasure, clustermatrix)
```

![fig. 1](figure/reproduced_figA1.png) 

The heatmap does not look exactly as [Fig. A.1](http://biostatistics.oxfordjournals.org/content/suppl/2006/04/21/kxj037.DC1/kxj037supp.pdf) in Johnson et al. This could be due to undocumented transformations preformed in Johnson et al., for example log-transformation or standardization. And there are different ways of calculating variation for a probe. The parameters used in the clustering could also be different. But the purpose of the figure was to show that a batch effect is present in the data. This is also achieved in the reproduced figure.
<BR/>  
<BR/>  
Next 3 data sets were made,
- **EB2**: Batches 1 and 2 adjusted with ComBat
- **EB3**: All the batches adjusted with ComBat
- **Batch3**: Only batch 3 (no adjustments)

```r
mat = datamatrix[,sampleannotation$Batch %in% c("1","2")]
EB2 = as.matrix( sva::ComBat(
            dat=mat, 
            batch=sampleannotation[colnames(mat),"Batch"], 
            mod=model.matrix( ~as.factor(sampleannotation[colnames(mat),"Type"])  ), 
            numCovs=NULL, 
            par.prior=useparprior, 
            prior.plots=FALSE))
rm(mat)
EB3 = as.matrix( sva::ComBat(
            dat=datamatrix, 
            batch=sampleannotation[colnames(datamatrix),"Batch"], 
            mod=model.matrix(~as.factor(sampleannotation[colnames(datamatrix),"Type"])), 
            numCovs=NULL, 
            par.prior=useparprior, 
            prior.plots=FALSE))
Batch3 = datamatrix[,sampleannotation$Batch=="3"]
```


Next we try to re-create Figure A.2

> A heatmap diagram of 770 genes from data set 2 after applying the EB batch adjustments.<cite> Johnson et al.

There is no description in how the 770 genes were selected so we use the same approach as in the reproduction of Figure A.1.

```r
variationmeasure = apply(EB3, 1, FUN= function(x){var(x)})
clustermatrix = EB3[order(variationmeasure, decreasing=TRUE),][1:770,]
aheatmap(clustermatrix, scale = "row", Rowv=T, Colv=T, 
         color = colorRampPalette(c("red", "black", "green"))(n = 299),
         annCol = sampleannotation[, c("Batch", "Type", "Cell")], 
         cellheight=NA, cellwidth=13,fontsize=12,   labRow=NA,
         main=paste("Reproduced Fig A2.",sep=""),border_color=NA,treeheight=c(0,100) )

rm(variationmeasure, clustermatrix)
```

![plot of chunk reproduced_figA2](figure/reproduced_figA2.png) 

Again the heatmap is not exactly as in Johnson et al, but the batch clustering is broken, and the samples cluster more by cell type and treatment type.

Now follows a few tests for differentially expressed probes.
> Differential expression was assessed using Welch’s t-test to determine the differential expression of RNAi versus control samples. EB2 produced at list of 86 significant genes at a false discovery (q-value) threshold of 0.05 (Storey and Tibshirani, 2003).<cite> Johnson et al.


```r
EB2_pvals = apply(EB2, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(EB2), "Type"]=="C"],
                                    x[sampleannotation[colnames(EB2), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(EB2_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 54660    15
```

Original number was **86**, reproduced number is **15**.


> The third batch alone produced a list of 37 significant genes using the same threshold. <cite> Johnson et al.


```r
Batch3_pvals = apply(Batch3, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(Batch3), "Type"]=="C"],
                                    x[sampleannotation[colnames(Batch3), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(Batch3_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 54657    18
```

Original number was **37**, reproduced number is **18**.


> Without any adjustment, combining these two batches produced a list of only 9 genes a q-value cutoff of 0.05<cite> Johnson et al.


```r
Batch12 = datamatrix[,sampleannotation$Batch %in% c("1","2")]
Batch12_pvals = apply(Batch12, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(Batch12), "Type"]=="C"],
                                    x[sampleannotation[colnames(Batch12), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(Batch12_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 54671     4
```

Original number was **9**, reproduced number is **4**.

> Welch’s t-test was also applied to EB3 to find differential expressed genes; yielding 1599 genes significant at a q-value cutoff of 0.05. <cite> Johnson et al.


```r
EB3_pvals = apply(EB3, 1 , 
                  FUN=function(x){t.test(
                                    x[sampleannotation[colnames(EB3), "Type"]=="C"],
                                    x[sampleannotation[colnames(EB3), "Type"]=="R"]
                                    )$p.value})
print(table(qvalue(EB3_pvals)$qvalue<0.05))
```

```
## 
## FALSE  TRUE 
## 53672  1003
```

Original number was **1599**, reproduced number is **1003**.

> Reducing the q-value threshold to 0.01 yielded 488 significant genes. <cite> Johnson et al.


```r
print(table(qvalue(EB3_pvals)$qvalue < 0.01))
```

```
## 
## FALSE  TRUE 
## 54330   345
```

Original number was **488**, reproduced number is **345**.

> ..decreasing the threshold further to 0.001 yielded 161 significant genes.<cite> Johnson et al.


```r
print(table(qvalue(EB3_pvals)$qvalue < 0.001))
```

```
## 
## FALSE  TRUE 
## 54572   103
```

Original number was **161**, reproduced number is **103**.

The reproduced numbers are not the same as the reported ones, but they are not completely off. This is likely due to undocumented steps performed by Johnson et al.

The above effort done in reproducing some of the results for data set 2 is done using non-log transformed data that also includes some negative values. In the next part of this report we will do an alternative analysis were batch is included in the analysis using limma instead of using ComBat. Since Limma needs log2 transformed data, the results will not be directly comparable to a t-test result using non-log values. Later we will also perform a few test using permutation (i.e many runs). For these reasons the rest of the analysis, including the significance test on ComBat adjusted data, will be performed on log2 transformed data with negative values set to 1 or removed (if in many samples for that probe). This is a fairly standard workflow and more practical for doing comparisons. In addition the parametric version of ComBat will be used.

First we make a data matrix with log2 values. And doing so the negative values will be floored and removed if in 

```r
log2data = datamatrix

# flooring to 1
log2data[log2data < 1] = 1

# take out data with to much low/missing values.
negativeprobesfilter = (rowSums(log2data > 1) >= (0.9 * ncol(log2data)))
log2data = log2data[negativeprobesfilter, ]

# quantilenormalize
log2data = normalizeBetweenArrays(log2(log2data))
```




```r

log2dataComBat = as.matrix(sva::ComBat(dat = log2data, batch = sampleannotation[colnames(log2data), 
    "Batch"], mod = model.matrix(~as.factor(sampleannotation[colnames(log2data), 
    "Type"])), numCovs = NULL, par.prior = useparprior, prior.plots = FALSE))
```

```
## Found 3 batches
## Found 1  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r

Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
fit = lmFit(log2dataComBat, design)
cont.matrix = makeContrasts(contrasts = "TypeR-TypeC", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
log2ComBat_res = eBayes(fit2)
rm(Type, design, fit, cont.matrix, fit2)
print(table(qvalue(log2ComBat_res$p.value[, 1])$qvalue < 0.05))
```

```
## 
## FALSE  TRUE 
## 51456   819
```

The number of genes below the significance threshold is somewhat lower than for the original the t-test, but that comparison is not relevant since the data is processed in different ways before the test. But for the explained practical purposes this is the result coming from the ComBat adjusted data.


### Analysis without ComBat
An alternative (and better?) way of handling batch effect is to include it in the model for the statistical test. This is possible in several tools and we choose to use the popular limma package (Smyth et al). 

We start with the data before ComBat adjustment. Limma works best with log transformed and between array normalized values, so first we set negative values to 1 and filter out probes that has this in more than half the samples. To ease later comparisons the same filter will also be applied to the EB3 data set. Then the test is run with the batch included as a blocking factor.

```r
Type = as.factor(sampleannotation$Type)
Block = as.factor(sampleannotation$Batch)
design = model.matrix(~0 + Type + Block)
fit = lmFit(log2data, design)
cont.matrix = makeContrasts(contrasts = "TypeR-TypeC", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
log2limmablock_res = eBayes(fit2)
rm(Type, Block, design, fit, cont.matrix, fit2)
print(table(qvalue(log2limmablock_res$p.value[, 1])$qvalue < 0.05))
```

```
## 
## FALSE  TRUE 
## 51878   397
```

Number of differentially expressed probes found with the alternative analysis is **397**.

Alternativly, limma´s <i>duplicateCorrelation</i> function can be used to treat batch as a random effect.

```r
Type = as.factor(sampleannotation$Type)
Block = as.factor(sampleannotation$Batch)
design = model.matrix(~0 + Type)
corfit <- duplicateCorrelation(log2data, design, block = Block)
fit <- lmFit(log2data, design, block = Block, correlation = corfit$consensus)
cont.matrix = makeContrasts(contrasts = "TypeR-TypeC", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
log2limmamixed_res = eBayes(fit2)
rm(Type, Block, design, fit, cont.matrix, fit2)
print(table(qvalue(log2limmamixed_res$p.value[, 1])$qvalue < 0.05, qvalue(log2limmablock_res$p.value[, 
    1])$qvalue < 0.05, dnn = c("duplicateCorrelation", "factor")))
```

```
##                     factor
## duplicateCorrelation FALSE  TRUE
##                FALSE 51848    60
##                TRUE     30   337
```

For this data set the two methods will give quite similar result. We will only use the first in subsequent comparisons.



### Consequences of ComBat use for the end result


Now we can compare the P-values for the two alternative procedures.

```r
hist(log2ComBat_res$p.value[, 1], border = "blue", main = "P-values, ComBat vs Limma", 
    breaks = 100, xlab = "p-value")
hist(log2limmablock_res$p.value[, 1], border = "red", add = T, breaks = 100)
legend("topright", legend = c("ComBat adjusted", "Limma adjusted"), text.col = c("blue", 
    "red"))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.svg) 

This plot shows that the p-values calculated from the ComBat adjusted data is skewed towards the low end compared to when batch is considered inside the statistical test in limma. However, the difference is not large and from this plot alone it could be argued that ComBat is just better.

Lets assume one aim of this study was to obtain a list of genes that differed between control and treatment with a significant cut-off of q<0.05. Then handling the batch effect with ComBat adjustment will yield a list of **819** whereas a list of **397** probes would have been closer to the truth. Also it is of interest to note that it is mostly the same probes that are found.

```r
table( qvalue(log2limmablock_res$p.value[,1])$qvalue<0.05 ,
         qvalue(log2ComBat_res$p.value[,1])$qvalue<0.05,
       dnn=c("limma", "ComBat"))
```

```
##        ComBat
## limma   FALSE  TRUE
##   FALSE 51455   423
##   TRUE      1   396
```

And if they were to use the top 1000 probes in a gene set test, they would have found many of the same genes regardless of handling batch effects by ComBat or limma.

```r
table( rank(log2limmablock_res$p.value[,1])  <=1000 & 
         rank(log2ComBat_res$p.value[,1])<=1000   )
```

```
## 
## FALSE  TRUE 
## 51412   863
```

Our conclusion is that for this study the error introduced by the use of ComBat would probably have a modest effect on the final result.

### Additional sanity checks

To substantiate that the result from the use of ComBat is less trustworthy than the alternative analysis we provide a few additional sanity checks.

First we use random numbers drawn from the same distribution regardless of batch or covariate but retaining the batch/covariate design.

```r
set.seed(100)
randomdata = log2data
randomdata[,] =rnorm(length(randomdata), mean=0, sd=1)
randomdataComBat = as.matrix( sva::ComBat(
            dat=randomdata, 
            batch=sampleannotation[colnames(randomdata),"Batch"], 
            mod=model.matrix(~as.factor(sampleannotation[colnames(randomdata),"Type"])), 
            numCovs=NULL, 
            par.prior=useparprior, 
            prior.plots=FALSE))
```

Limma is then used in 3 ways
- On the ComBat adjusted random numbers
- On the random numbers including batch as a blocking factor
- On the random numbers ignoring batch information

```r

Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  

fit = lmFit(randomdataComBat, design)
fit2 = contrasts.fit(fit, cont.matrix)
randomdataComBat_pvalues = eBayes(fit2)$p.value[,1]


fit = lmFit(randomdata, design)
fit2 = contrasts.fit(fit, cont.matrix)
randomdata_pvalues = eBayes(fit2)$p.value[,1]



Block = as.factor(sampleannotation$Batch)
design = model.matrix(~0+Type+Block)
fit = lmFit(randomdata, design)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design) 
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
randomdatalimmablock_pvalues = fit2$p.value[,1]

```

Plotting the p-values

```r
hist(randomdataComBat_pvalues, border = "blue", main = "P-values, Random numbers", 
    breaks = 100, xlab = "p-value")
hist(randomdata_pvalues, border = "black", add = T, breaks = 100)
hist(randomdatalimmablock_pvalues, border = "red", add = T, breaks = 100)
legend("topright", legend = c("ComBat adjusted", "Limma adjusted", "No adjustment"), 
    text.col = c("blue", "black", "red"))
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.svg) 

The p-values for ComBat adjusted random numbers are slightly enriched for low p-values. Indicating that some of the enrichment seen in the real data is also false.
<br/>  
<br/>  
Another sanity check is to subset the real data into fictive batches. The idea behind this check stems from the observation that the EB2 data set contains the same number of samples with the same batch/covariate balance as the third batch alone.  

> The third batch was used for comparison against the EB2 analysis results because it was an identical experiment to EB2 other than the fact that it was conducted in a single batch. <cite> Johnson et al.

Thus it is interesting to inspect the p-value from these two identical experiments, EB2 (batch 1 and 2 adjusted with ComBat) and batch 3 (no ComBat adjustment). These p-values are already computed above. 

```r
# need to subset the log2data and do the ComBat adjustment and do the limma
# tests

# subset log2data
mat = log2data[, sampleannotation$Batch %in% c("1", "2")]
log2EB2 = as.matrix(sva::ComBat(dat = mat, batch = sampleannotation[colnames(mat), 
    "Batch"], mod = model.matrix(~as.factor(sampleannotation[colnames(mat), 
    "Type"])), numCovs = NULL, par.prior = useparprior, prior.plots = FALSE))
```

```
## Found 2 batches
## Found 1  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r
rm(mat)
log2Batch3 = log2data[, sampleannotation$Batch == "3"]

# run limma
Type = as.factor(sampleannotation[colnames(log2EB2), "Type"])
design = model.matrix(~0 + Type)
fit = lmFit(log2EB2, design)
cont.matrix = makeContrasts(contrasts = "TypeR-TypeC", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
log2EB2ComBat_pvals = eBayes(fit2)$p.value
rm(Type, design, fit, cont.matrix, fit2)

Type = as.factor(sampleannotation[colnames(log2Batch3), "Type"])
design = model.matrix(~0 + Type)
fit = lmFit(log2Batch3, design)
cont.matrix = makeContrasts(contrasts = "TypeR-TypeC", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
log2Batch3ComBat_pvals = eBayes(fit2)$p.value
rm(Type, design, fit, cont.matrix, fit2)

# plot the two p-value distibutions
hist(log2EB2ComBat_pvals, border = "blue", main = "P-values, real data", breaks = 100, 
    xlab = "p-value")
hist(log2Batch3ComBat_pvals, border = "red", add = T, breaks = 100)
legend("topright", legend = c("ComBat adjusted batch 1 and 2", "Batch 3"), text.col = c("blue", 
    "red"))
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.svg) 

The 2-batch experiment is able to retrieve more genes as significant(21 for q<0.05) than running all in one batch(29 for q<0.05) for the same sample size. This is also observed by Johnson et al., but not commented. We claim that the apparent better result for the two-batch experiment is a consequence of the use of ComBat. A sanity check supporting this claim is to split batch 3 into two imaginary batches with the same design as for batch 1 and 2, and then look at the p-value distribution compared with the real batch 3. In order to obtain some robustness this test in performed for 10 permutations of different ways of sub-setting batch 3 into 2 fictive batches. 



```r

# print the batch/covariate balance for the real batch 1 and batch 2
print(table(sampleannotation[sampleannotation$Batch %in% c("1", "2"), c("Batch", "Type")]))

fictivebatches_pvalcounts=list()
runs = 10
for(i in 1:runs)
{
permdata = log2Batch3
permdata_annot = sampleannotation[dimnames(permdata)[[2]],]
permdata_annot$Batch="4"
permdata_annot$Batch[permdata_annot$Type=="C"][sample(1:sum(permdata_annot$Type=="C"), 4)] = "5"
permdata_annot$Batch[permdata_annot$Type=="R"][sample(1:sum(permdata_annot$Type=="R"), 3)] = "5"

# This permutations balance
print(table(permdata_annot[, c("Batch", "Type")]))

permdataComBat = as.matrix( sva::ComBat(
            dat=permdata, 
            batch=permdata_annot[,"Batch"], 
            mod=model.matrix( ~as.factor(permdata_annot[,"Type"])  ), 
            numCovs=NULL, 
            par.prior=TRUE, 
            prior.plots=FALSE))


Type = as.factor(permdata_annot[, "Type"])
design = model.matrix(~0 + Type)
fit = lmFit(permdataComBat, design)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
pvals = eBayes(fit2)$p.value

fictivebatches_pvalcounts[[i]] = hist(pvals, plot=FALSE, breaks=100)$counts
}
```



The p-values for the 10 splits of Batch 3 adjusted by ComBat is compared to the p-values for the not-split(and not ComBat adjusted) Batch 3.

```r
plot((1:20)/100, hist(log2Batch3ComBat_pvals, plot=FALSE, breaks=100)$counts[1:20],
     col="red",  main="P-values, fictive batches", 
      xlab="p-value", ylab="Frequency", type="l", lwd=2, 
      ylim=c(0, max(unlist(fictivebatches_pvalcounts))))
for(i in 1:length(fictivebatches_pvalcounts))
{
  lines((1:20)/100,fictivebatches_pvalcounts[[i]][1:20], col="blue")
}
legend("topright", legend=c("ComBat adjusted fictive-batches of batch 3", "Batch 3 not adjusted"), text.col=c("blue", "red"))
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.svg) 

The two-batch approach seem to consistently outperform the all-in-one-batch version. 




### References


Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127.

Storey, J. D. and Tibshirani, R. (2003) Proc Natl Acad Sci U S A, 100, 9440-5.

Raivo Kolde (2013). pheatmap: Pretty Heatmaps. R package version 0.7.7. http://CRAN.R-project.org/package=pheatmap

Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD.(2012) The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics. 2012 Mar 15;28(6):882-3.

Smyth, GK (2005). Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions
  using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages
  397-420.
  
  R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing,
  Vienna, Austria. URL http://www.R-project.org/

  Yihui Xie (2013). knitr: A general-purpose package for dynamic report generation in R. R package version 1.5.

  Yihui Xie (2013) Dynamic Documents with R and knitr. Chapman and Hall/CRC. ISBN 978-1482203530

  Yihui Xie (2013) knitr: A Comprehensive Tool for Reproducible Research in R. In Victoria Stodden, Friedrich Leisch and
  Roger D. Peng, editors, Implementing Reproducible Computational Research. Chapman and Hall/CRC. ISBN 978-1466561595
  
  RStudio Team (2012). RStudio: Integrated Development for R. RStudio, Inc., Boston, MA URL http://www.rstudio.com/.



```r
sessionInfo()
```

```
R version 3.0.2 (2013-09-25)
Platform: x86_64-apple-darwin10.8.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] statmod_1.4.19     NMF_0.20.5         Biobase_2.22.0    
 [4] BiocGenerics_0.8.0 cluster_1.15.2     rngtools_1.2.4    
 [7] pkgmaker_0.20      registry_0.2       limma_3.18.13     
[10] qvalue_1.36.0      sva_3.8.0          mgcv_1.7-29       
[13] nlme_3.1-117       corpcor_1.6.6      pheatmap_0.7.7    
[16] knitr_1.5         

loaded via a namespace (and not attached):
 [1] codetools_0.2-8    colorspace_1.2-4   digest_0.6.4      
 [4] doParallel_1.0.8   evaluate_0.5.5     foreach_1.4.2     
 [7] formatR_0.10       ggplot2_0.9.3.1    grid_3.0.2        
[10] gridBase_0.4-7     gtable_0.1.2       iterators_1.0.7   
[13] lattice_0.20-29    MASS_7.3-33        Matrix_1.1-3      
[16] munsell_0.4.2      plyr_1.8.1         proto_0.3-10      
[19] RColorBrewer_1.0-5 Rcpp_0.11.1        reshape2_1.4      
[22] scales_0.2.4       stringr_0.6.2      tcltk_3.0.2       
[25] tools_3.0.2        xtable_1.7-3      
```


generation ended 2014-06-11 14:21:26. Time spent 8 minutes .

