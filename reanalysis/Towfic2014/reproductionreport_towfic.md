Partial reproduction and alternative analysis of "Comparing the Biological Impact of Glatiramer Acetate with the Biological Impact of a Generic"
========================================================

2014-05-23 23:19:35

### Overview
This report aims to show that the use of the statistical tool ComBat from [Johnson et al.](http://biostatistics.oxfordjournals.org/content/8/1/118.abstract) led to false results in the analysis performed in [Towfic et al.'s ](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0083757) "Comparing the Biological Impact of Glatiramer Acetate with the Biological Impact of a Generic".

The generation of this report is a side product of  "Methods that remove batch effects while retaining group differences may lead
to exaggerated confidence in downstream analyses", V. Nygaard, E. A. Rødland, E. Hovig, manuscript in preparation. It consist of some of the analyses performed working with the manuscript that were to detailed to be included. Nevertheless it should be useful for the specially interested. This report is not peer-reviewed.

The format of this report is html with text, r-code and plots intermingled and made from a rmarkdown file in R-studio with knitr. It is meant to be reproducible if run in conjunction with a few accompanying files from the github repository.

This document has five main parts
- Getting and formatting the data
- Reproduce some of the results to show that we are working on the same data and analysis workflow
- Remove the use of ComBat and perform a similar analysis with an alternative established tool
- Estimate the error introduced when ComBat is used and the consequences for the conclusion of the study
- Perform a few more sanity checks to substantiate that the difference in results for the two above analysis is mainly introduced by ComBat






```r
includelibs = c("Biobase", "GEOquery", "sva", "limma", "preprocessCore", "statmod")
lapply(includelibs, require, character.only = T)
```

```
## Warning: package 'statmod' was built under R version 3.0.3
```

```r
source("../../commonscripts/helperfunctions.r")
source("helperfunctions_towfic.r")
```

The above libraries and helper script files are needed in order to reproduce this report.

### Getting the data

The important sample information is described in Table S1 from Towfic et al. and its usage in ComBat is described briefly in the methods;
>  Each microarrays chip designation was supplied a batch label; there were 18 batches in all. The labels for the treatments (i.e. drug product, reference standard) were added as covariates.<cite> Towfic et al.

Table S1 does have a "Chip"- column, unfortunately there is no dedicated "treatment"-column.
Communication with the corresponding author yielded this explanation
> ..and the only thing we used as a covariate was the unique treatment names themselves (e.g. "Medium"" or "RS").

Based on these descriptions and the annotation for the GEO deposit, we compiled a more complete sample annotation file (sampleannotation.csv).


```r
sampleannotation = read.table("data/sampleannotation.csv", sep = "\t", header = TRUE, 
    stringsAsFactors = FALSE)
sampleannotation$code = make.names(sampleannotation$code)
sampleannotation$chip = as.character(sampleannotation$chip)
dimnames(sampleannotation)[[1]] = sampleannotation$code
head(sampleannotation)
```

```
##       number  code covariate       chip slot geogroup geoaccession
## RS.16      1 RS.16        RS 4634633002    A Verified    GSM996691
## DP.34      2 DP.34        DP 4634633002    B Verified    GSM996731
## DP.19      3 DP.19        DP 4634633002    C Verified    GSM996716
## DP.32      4 DP.32        DP 4634633002    D Verified    GSM996729
## DP.04      5 DP.04        DP 4634633002    E Verified    GSM996701
## MBP        6   MBP       MBP 4634633002    F    other    GSM996658
```

The "covariate"-column is made based on the "code"-column, but might not match 100% to what were actually used for all samples. The last two columns are from the GEO deposit and reveal the samples that do not have data, which we will remove from the data shortly.

The naming convention seem to differ slightly between the Table S1 and the text in Towfic et al. This is our interpretation of the main covariate labels and its corresponding name in the text.
- **DP** referred to as GA (but not GA as in table S1)
- **N** referred to as "generic"
- **M** referred to as "medium"
- **RS** referred to as "reference standard""


```r
# take out 3 samples that are not assign to a geoaccession. Failed QC?
sampleannotation = sampleannotation[!is.na(sampleannotation$geoaccession), ]
```

The batch/covariate design shows many batches and covariate groups.

```r
table(sampleannotation[, c("chip", "covariate")])
```

```
##             covariate
## chip         A C DM DP E GA M MANITOL MBP N RS S TV
##   4634633002 0 0  0  4 0  0 0       0   1 0  1 0  0
##   4634633043 0 0  1  0 0  1 1       0   0 2  1 0  0
##   4634633053 0 0  1  1 0  0 0       0   0 1  1 0  0
##   4637105025 0 0  0  5 0  0 0       1   0 0  0 0  0
##   4682416019 0 0  0  4 0  0 0       0   0 0  2 0  0
##   4682416043 0 0  0  1 0  2 0       0   0 1  2 0  0
##   4682416087 0 0  1  0 0  0 1       0   0 0  1 0  3
##   4682416090 0 0  0  4 0  0 1       0   0 1  0 0  0
##   4763128059 0 0  1  2 0  0 1       0   0 1  1 0  0
##   4763128060 0 0  1  2 0  0 1       0   0 1  1 0  0
##   4763646004 0 0  0  1 0  0 1       0   0 0  4 0  0
##   4763646012 0 1  0  1 0  0 1       0   0 2  1 0  0
##   4763646026 0 1  0  1 0  0 1       0   0 2  1 0  0
##   4763646030 0 0  0  2 0  0 2       0   0 0  2 0  0
##   5216898004 1 0  0  2 0  0 1       0   0 0  1 0  0
##   5216898008 0 0  0  2 0  0 1       0   0 0  1 1  1
##   5216898012 2 0  0  1 1  0 1       0   0 0  1 0  0
##   5216898024 0 0  0  1 1  0 1       0   0 0  1 1  1
```

A look at the primary comparison, "DP"(GA) and "N"(generic) reveals a lack of balance.

```r
table(sampleannotation[sampleannotation$covariate %in% c("DP", "N"), c("chip", 
    "covariate")])
```

```
##             covariate
## chip         DP N
##   4634633002  4 0
##   4634633043  0 2
##   4634633053  1 1
##   4637105025  5 0
##   4682416019  4 0
##   4682416043  1 1
##   4682416090  4 1
##   4763128059  2 1
##   4763128060  2 1
##   4763646004  1 0
##   4763646012  1 2
##   4763646026  1 2
##   4763646030  2 0
##   5216898004  2 0
##   5216898008  2 0
##   5216898012  1 0
##   5216898024  1 0
```


> The microarray data have been deposited in the Gene Expression Omnibus, under accession number [GSE40566](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40566)

The processed data from the above GEO deposit ("Series Matrix File") is from another publication which did not processes the data as described in Towfic et al. (personal communication). The lesser processed data in GEO linked to as [GSE40566_non_normalized.txt.gz](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz) in the "Supplementary file" section was thus used in this report.

Probe annotation was also found in the "Supplementary file" section as [GSE40566_RAW.tar](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file).

A recap of the important input files used in this report:
- **sampleannotation.csv** holds the important sample/batch/covariate assignment. Compiled based on information in Towfic et al, the geo deposit and personal communication.
- **GSE40566_non_normalized.txt** the measurements in a sample vs probe matrix. From the [GEO deposit](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz). Not in github.
- **GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt** the probe annotation needed for the data row to probe to gene matching. From the GEO deposit inside the [GSE40566_RAW.tar](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file). Not in github.

Reading the data.

```r
if (downloaddata) {
    temp = tempfile()
    download.file(url = "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz", 
        destfile = temp, mode = "wb")
    rawdata = read.table(temp, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    unlink(temp)
    
    temp = tempfile()
    download.file(url = "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file", 
        destfile = temp, mode = "wb")
    tardirtemp = tempdir()
    untar(temp, exdir = tardirtemp)
    rawannotation = read.table(paste(tardirtemp, "/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt.gz", 
        sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, skip = 8, 
        comment.char = "", quote = "", fill = TRUE)
    
    unlink(temp)
    unlink(tardirtemp, recursive = TRUE)
} else {
    
    # if download did not work change downloaddata to FALSE download data from
    # the GEO deposit unpack and place files in a folder named 'not_in_github'
    rawdata = read.table("not_in_github/GSE40566_non_normalized.txt", sep = "\t", 
        header = TRUE, stringsAsFactors = FALSE)
    
    # the probe annotation file found inside GSE40566_RAW.tar
    rawannotation = read.table("not_in_github/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt", 
        sep = "\t", header = TRUE, stringsAsFactors = FALSE, skip = 8, comment.char = "", 
        quote = "", fill = TRUE)
}
```


Some formatting is needed to match the probe information to the measurements and the measurements to the sample information.

```r
# annotation file has two separate tables, one for experimental and one at the end for controls.
# splitting and fixing the tables:
tmp =  rawannotation$Species=="Mus musculus"
experimentalannot = rawannotation[tmp,]
experimentalannot$Array_Address_Id = as.numeric(experimentalannot$Array_Address_Id)
controlannot = rawannotation[!tmp,]
dimnames(controlannot)[[2]] = rawannotation[rawannotation[,2]=="Array_Address_Id",]
controlannot$Array_Address_Id = suppressWarnings(as.numeric(controlannot$Array_Address_Id))
controlannot = controlannot[!is.na(controlannot$Array_Address_Id),]
controlannot=controlannot[,1:6]
probeannotation = merge(experimentalannot, controlannot, all=TRUE )
#dim(probeannotation)
rm(tmp, experimentalannot, controlannot)

probeannotation = probeannotation[!duplicated(probeannotation$Array_Address_Id),]
probeannotation = probeannotation[probeannotation$Array_Address_Id %in% rawdata$ID_REF, ] # 
probeannotation$Symbol=tolower(probeannotation$Symbol)
dimnames(probeannotation)[[1]] = probeannotation$Probe_Id
#dim(probeannotation)

#sort and filter probe and data similar.
datamatrix_raw = as.matrix(rawdata[,-1])
datamatrix_raw = datamatrix_raw[match( probeannotation$Array_Address_Id , rawdata$ID_REF), ]
dimnames(datamatrix_raw)[[1]] = probeannotation$Probe_Id
#dim(datamatrix_raw)
#dim(probeannotation)

#and match data to samples.
#table(sampleannotation$code %in% dimnames(datamatrix_raw)[[2]])# check
#table(dimnames(datamatrix_raw)[[2]] %in% sampleannotation$code)# check
datamatrix_raw = datamatrix_raw[, match(sampleannotation$code , dimnames(datamatrix_raw)[[2]])]
```


Several tables with results are presented in the Supporting Information in Towfic et al. We are aiming to reproduce two of these. Unfortunately they were only presented in a pdf-format not easily parsed by a computer. We had to resort to a adhoc method of cut-and-paste from pdf into text files which were somewhat polluted by pdf-formatting code. For some tables about 10% of the rows will be lost, but the rest will suffice for our purpose.


```r
# Table S2
table_s2 = read.table("data/table_s2.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    fill = TRUE)
table_s2 = as.matrix(table_s2)
table_s2[, 3:6] = as.numeric(table_s2[, 3:6])
```

```
## Warning: NAs introduced by coercion
```

```r
# not able to paste the pdf whitout a lot of gibberish clutter the data.
# Some probes are lost!
a = rowSums(is.na(table_s2[, 3:6])) > 0
print(paste("Lost rows table s2: ", sum(a)))
```

```
## [1] "Lost rows table s2:  91"
```

```r
table_s2 = data.frame(table_s2[!a, ], stringsAsFactors = FALSE)
# data.frame makes the columns characters again.
for (n in 3:6) {
    table_s2[, n] = as.numeric(table_s2[, n])
}

# Table S5
table_s5 = read.table("data/table_s5.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    fill = TRUE)
tmp = as.matrix(table_s5)
tmp[, 3:14] = as.numeric(tmp[, 3:14])
```

```
## Warning: NAs introduced by coercion
```

```r
# not able to paste the pdf whitout a lot of gibberish clutter the data.
# Some probes are lost!
a = rowSums(is.na(tmp[, 3:14])) > 0
print(paste("Lost rows table s5: ", sum(a)))
```

```
## [1] "Lost rows table s5:  15"
```

```r
table_s5 = data.frame(table_s5[!a, ], stringsAsFactors = FALSE)
table_s5$Fold_Change = as.numeric(table_s5$Fold_Change)
```



### Reproduce some of the original results

Towfic et al. performed many different test on the data and it is outside the scope of this report to reproduce all of their results. We focus on the key part of testing for differentially expressed genes between "GA" (DP) and "generic"(N) as described in Table S2 and Table S5. But before those tests we have to preprocess according to the description.

> Starting with background-corrected bead-level signals, we quantile normalized the extracted data for all samples across all 46,547 probes via the preprocessCore package in R.


```r
datamatrices = list()
datamatrices[["real_raw"]] = datamatrix_raw
if (debug) datamatrices[["real_raw"]] = datamatrix_raw[1:1000, ]

datamatrices[["real_qnorm"]] = normalize.quantiles(datamatrices[["real_raw"]])
# used in paper, but seems to do the same as the limma version. Except it
# loses dimnames
dimnames(datamatrices[["real_qnorm"]]) = dimnames(datamatrices[["real_raw"]])
```


> We then corrected for batch variation with ComBat [17] as implemented in the SVA package of R [18]. Each microarrays chip designation was supplied a batch label; there were 18 batches in all. The labels for the treatments (i.e. drug product, reference standard) were added as covariates. 


```r
combatmod = model.matrix(~as.factor(sampleannotation$covariate))
datamatrices[["real_combat_covariates"]] = as.matrix(ComBat(dat = datamatrices[["real_qnorm"]], 
    batch = sampleannotation$chip, mod = combatmod, numCovs = NULL, par.prior = TRUE, 
    prior.plots = FALSE))
```

```
## Found 18 batches
## Found 12  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

**Table S2**

First we will try to reproduce Table S2.

> Genes utilized for the tolerance method illustrated in Figure 1B. (figure text)

> This standard of comparison was constructed by first identifying the top 1000 probes by absolute fold change of reference standard compared to the medium (Table S2). The list includes both upregulated and downregulated probes compared to medium. Probes were filtered such that ones upregulated by reference standard needed to have an average reference standard expression of 6.00 or higher and ones downregulated by reference standard needed to have an average medium expression of 6.00 or higher. 

This is how the table looks

```r
head(table_s2)
```

```
##             ID         Gene AVG_Medium AVG_Reference_Standard AVG_GA
## 1 ILMN_2685712         IFNG      7.152                 11.255 11.163
## 2 ILMN_1247309          IL3      5.823                  9.406  9.384
## 3 ILMN_2791459         IFNG      6.636                  9.829  9.711
## 4 ILMN_1215862        CXCL9      7.342                  9.817  9.757
## 5 ILMN_2595732 LOC100046232      9.249                 11.509 11.547
## 6 ILMN_2931334          IL4      6.172                  8.367  8.425
##   AVG_generic
## 1      10.668
## 2       9.422
## 3       9.304
## 4       9.525
## 5      11.767
## 6       8.496
```



```r
thiscovariatename = "covariate"
thisdata = datamatrices[["real_combat_covariates"]]

# the mean uses log2 numbers and will be a geometric mean.
meanM = rowMeans(thisdata[, sampleannotation[, "covariate"] == "M"])
meanRS = rowMeans(thisdata[, sampleannotation[, "covariate"] == "RS"])
meanDP = rowMeans(thisdata[, sampleannotation[, "covariate"] == "DP"])
meanN = rowMeans(thisdata[, sampleannotation[, "covariate"] == "N"])
rm(thisdata)
foldchange_RS_vs_M = meanRS - meanM
foldchange_RS_vs_M = foldchange_RS_vs_M[meanM > 6 | meanRS > 6]
top1000_RS_vs_M = names(foldchange_RS_vs_M[order(abs(foldchange_RS_vs_M), decreasing = TRUE)])[1:1000]
table(table_s2$ID %in% top1000_RS_vs_M)
```

```
## 
## FALSE  TRUE 
##    58   851
```

Out of Table S2's **909** probes, **851** were reproduced. Scatter plots of the actual values shows a correlation.


```r
# scatter plot correlation of the 4 columns in table_S2 vs. the reproduced numbers.
par(mfrow=c(2, 2))
plot(table_s2$AVG_Medium, meanM[table_s2$ID], pch=".")
plot(table_s2$AVG_Reference_Standard, meanRS[table_s2$ID], pch=".")
plot(table_s2$AVG_GA, meanDP[table_s2$ID], pch=".")
plot(table_s2$AVG_generic, meanN[table_s2$ID], pch=".")
```

![plot of chunk reproduced_tables2](figure/reproduced_tables2.png) 

Original values from table S2 on X-axis, reporduced values on Y-axis. The values in Table S2 is almost reproduced.


**Table S5**

This is how the table looks. It is transposed and shows only the top 3 genes for visualization purposes.

```r
t(table_s5[1:3, ])
```

```
##                   1              2              3             
## Probe             "ILMN_2685712" "ILMN_2595857" "ILMN_2791459"
## Gene              "IFNG"         "ANKRD37"      "IFNG"        
## Fold_Change       "1.41"         "1.36"         "1.33"        
## t_test_p          "0.0002995"    "0.0002995"    "0.0002995"   
## t_test_adjp       "0.01489"      "0.01489"      "0.01489"     
## snr_p             "0.0002995"    "0.0002995"    "0.0002995"   
## snr_adjp          "0.01335"      "0.01335"      "0.01335"     
## limma_p           "2.36e-07"     "5.82e-12"     "1.99e-06"    
## limma_adjp        "7.180e-05"    "1.690e-08"    "3.756e-04"   
## anova_p           "4.01e-07"     "1.23e-08"     "6.42e-06"    
## anova_adjp        "1.39e-04"     "1.19e-05"     "1.08e-03"    
## signif_parametric "1"            "1"            "1"           
## wilcoxon_p        "2.065e-05"    "6.035e-07"    "1.972e-04"   
## wilcoxon_adjp     "0.0043099"    "0.0005733"    "0.0154771"
```


> Comparison of expression in GA to expression in GA for each probe, including fold change, ANOVA, LIMMA with background subtraction, comparative marker selection by signal-to-noise ratio, comparative marker selection by t-test, and the Wilcoxon non-parametric method. <cite> Figure text, Towfic et al.

There must be a typo in the figure text, they meant "Comparison of expression in GA to expression in **generic** for each probe, ...".

> To find differentially expressed probes between generic and GA, we utilized various statistical tests at the probe level and merged the results across the different methods. <cite> Towfic et al.

How the fold change was calculated is not explained. A common way is to take the difference of geometric means that were calculated for Table S2 above.

```r
foldchange_DP_N = meanDP - meanN
```


These reproduced fold changes can be compared for the probes listed in table S5.


```r
par(mfrow=c(1, 2))
plot(table_s5$Fold_Change, foldchange_DP_N[table_s5$Probe],
     pch=".", xlab="Original Table S5", ylab="Reproduced",
     main="FC comparison")

thiscolors = c("black", "red")
hist(foldchange_DP_N, breaks=100, freq=F, border=thiscolors[1], main="Reproduced FC for DP vs. N")
hist(foldchange_DP_N[table_s5$Probe], breaks=100, add=T, freq=F, border=thiscolors[2])
legend("topright", legend=c("All probes", "Table S5 probes"), text.col=thiscolors)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

The first plot reveals that the reproduced fold changes are much more compressed than the ones in Table S5.
There seems to two linear relationship, one for negative FC's and another for positive FC's. They seems to have different intercepts but comparable slope. The second plot shows that the reproduced fold changes for the probes listed in Table S5 are clearly more extreme than the rest.  

The fold changes in table 5S do not match the average expression values listed in table S2 either as a the same plot will show;


```r
# table S2 has the AVG for both DP and N listed
fc_table_s2 = table_s2$AVG_GA - table_s2$AVG_generic
names(fc_table_s2) = table_s2$ID
x = table_s5$Fold_Change
y = fc_table_s2[table_s5$Probe]
lim =c( min(   c(min(x, na.rm=T), min(y, na.rm=T))  ), max( c(max(x, na.rm=T), max(y, na.rm=T)) ))
plot(x, y , ylim=lim, xlim=lim, pch="." ,  xlab="Table_S5", ylab="Table_S2",
     main="FC comparison, given values.") #  almost perfect match with offset.
rm(x, y, lim)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 

There seems to be some inconsistency between between the fold changes in table S5 and the reproduced ones, however it is the confidence i.e p-values that are affected by the use of ComBat. 

The rest of table S5 is p-values and we will try to reproduce the ones from limma.

> Next, we utilized Linear Models for Microarray (LIMMA) data analysis [23], [24] R package, part of the Bioconductor framework [25], to compare generic and GA samples, fitting a linear model that adjusts for fixed effect from medium (Effect = (GA-generic)  (generic-Medium)). The coefficients for the linear model were tested for significance using a modified t-test (taking into account standard deviation)...<cite> Towfic et al.

The formula "Effect = (GA-generic)  (generic-Medium)" seems strange.


```r
group = factor(sampleannotation$covariate)
design = model.matrix(~0 + group)
fit = lmFit(datamatrices[["real_combat_covariates"]], design)
cont.matrix = makeContrasts ( contrasts="(groupDP-groupN) - (groupN-groupM)", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
limma_ret_org = eBayes(fit2)
limma_p_org = limma_ret_org$p.value
names(limma_p_org) = names(limma_ret_org$Amean)
par(mfrow=c(1, 2))
plot(table_s5[,"limma_p"], limma_p_org[table_s5$Probe], xlab="Original Table S5", 
     ylab="Reproduced", main="LIMMA P-values from table S5 vs. reproduced")
thiscolors = c("black", "red")
hist(limma_p_org, breaks=100, freq=T, border=thiscolors[1], main="Reproduced LIMMA p-values for DP vs. N")
hist(limma_p_org[table_s5$Probe], breaks=100, add=T, freq=T, border=thiscolors[2])
legend("topright", legend=c("All probes", "Table S5 probes"), text.col=thiscolors)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 

A lot of the reproduced p-values are high for the probes in Table S5. A more straightforward contrast of Effect = GA-generic give more similar p-values as in Table S5, and for the rest of this report this is what will be used.


```r
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
limma_ret_alt = eBayes(fit2)
limma_p_alt = limma_ret_alt$p.value
names(limma_p_alt) = names(limma_ret_alt$Amean)
par(mfrow=c(1, 2))
plot(table_s5[,"limma_p"], limma_p_alt[table_s5$Probe], xlab="Original Table S5", 
     ylab="Reproduced", main="LIMMA P-values from table S5 vs. reproduced")
thiscolors = c("black", "red")
hist(limma_p_alt, breaks=100, freq=T, border=thiscolors[1], main="Reproduced LIMMA p-values for DP vs. N")
hist(limma_p_alt[table_s5$Probe], breaks=100, add=T, freq=T, border=thiscolors[2])
legend("topright", legend=c("All probes", "Table S5 probes"), text.col=thiscolors)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 

The reproduced p-values for the Table S5 probes are not exactly as the listed ones, but still mostly low. In addition to the limma test, Towfic et al. performed similar tests for significant probes between generic and GA using several other methods. We tried several of them as well, but have skipped them in this report since they showed a similar pattern as for the limma p-values.

This concludes the reproduction. We managed to reproduce the average expressions listed in Table S2 pretty well. The fold changes and p-values in Table S5 were with some reasonable interpretation reproduced. The results can now be compared to the same analysis using an alternative to ComBat.



### Analysis without ComBat and consequences on the results

Table S2 has average expression values for the conditions, but includes no values about how much trust we can have in those average expressions (i.e significance tests, confidence intervals etc.).

Table S5 on the other hand consists of results from significance tests that are strongly influenced by ComBat's data transformation. An alternative way of handling batch effect compared to ComBat's is to include it in the statistical test. This can be done in limma


```r
group = factor(sampleannotation$covariate)
block = factor(sampleannotation$chip)
design = model.matrix(~0+group+block)
fit = lmFit(datamatrices[["real_qnorm"]], design)
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
#limma_ret_woc = eBayes(fit2)
#limma_p_woc = limma_ret_woc$p.value
#names(limma_p_woc) = names(limma_ret_woc$Amean)
limma_p_woc =  eBayes(fit2)$p.value
rm(design, group, block, fit, cont.matrix, fit2)
table(p.adjust(limma_p_alt , method="fdr")<0.05, p.adjust(limma_p_woc , method="fdr")<0.05, dnn=c("ComBat adj", "LIMMA adj"))
```

```
##           LIMMA adj
## ComBat adj FALSE  TRUE
##      FALSE 44493     0
##      TRUE   1735     9
```

The number of significant probes using the given FDR threshold of 0.05 when batch effect is handled by LIMMA is 
**9**. This is much lower compared to when ComBat is applied, 
**1744**.  The p-value distribution is still skewed, but not as much.

Alternatively, limma´s <i>duplicateCorrelation</i> function can be used to treat batch as a random effect. This could arguably be a more powerful approach, and a short investigation is warranted.

```r

group = factor(sampleannotation$covariate)
block = factor(sampleannotation$chip)
design = model.matrix(~0 + group)

corfit <- duplicateCorrelation(datamatrices[["real_qnorm"]], design, block = block)
fit <- lmFit(datamatrices[["real_qnorm"]], design, block = block, correlation = corfit$consensus)
cont.matrix = makeContrasts(contrasts = "groupDP-groupN", levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
# limma_ret_woc = eBayes(fit2)
limma_p_mixed = eBayes(fit2)$p.value
```

```r
thiscolors = c("red", "blue", "black")
hist(limma_p_alt, breaks=100, freq=T, border=thiscolors[1], main="Alternative P-value distributions for DP vs. N", xlab="P-values")
hist(limma_p_woc, breaks=100, add=T, freq=T, border=thiscolors[2])
hist(limma_p_mixed, breaks=100, add=T, freq=T, border=thiscolors[3])
legend("topright", legend=c("ComBat adjusted", "LIMMA adjusted as factor", "LIMMA adjusted as random effect"), text.col=thiscolors)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 


For this data set the two LIMMA methods will give quite similar results, at least compared to the ComBat adjusted approach. We will only use the first, "batch as factor", in subsequent comparisons.


### Additional sanity checks

To substantiate that the result from the use of ComBat is less trustworthy than the alternative analysis we provide a couple of additional sanity checks.

First we use random numbers drawn from the same distribution(mean=0, sd=1) regardless of batch or covariate but retaining the batch/covariate design.

```r
set.seed(100)
datamatrices[["random_raw"]] = matrix(rnorm(length(datamatrices[["real_raw"]]), 
    mean = 0, sd = 1), nrow = nrow(datamatrices[["real_raw"]]), ncol = ncol(datamatrices[["real_raw"]]), 
    dimnames = dimnames(datamatrices[["real_raw"]]))
datamatrices[["random_combat_covariates"]] = as.matrix(ComBat(dat = datamatrices[["random_raw"]], 
    batch = sampleannotation$chip, mod = combatmod, numCovs = NULL, par.prior = TRUE, 
    prior.plots = FALSE))
```

```
## Found 18 batches
## Found 12  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```


limma is then used in 3 ways
- On the ComBat adjusted random numbers
- On the random numbers ignoring batch information
- On the random numbers including batch as a blocking factor



```r
group = factor(sampleannotation$covariate)
design = model.matrix(~0 + group)
cont.matrix = makeContrasts(contrasts = "groupDP-groupN", levels = design)


fit = lmFit(datamatrices[["random_combat_covariates"]], design)
fit2 = contrasts.fit(fit, cont.matrix)
limma_p_rand_combat = eBayes(fit2)$p.value[, 1]

fit = lmFit(datamatrices[["random_raw"]], design)
fit2 = contrasts.fit(fit, cont.matrix)
limma_p_rand_nocombat = eBayes(fit2)$p.value[, 1]


block = as.factor(sampleannotation$chip)
design = model.matrix(~0 + group + block)
cont.matrix = makeContrasts(contrasts = "groupDP-groupN", levels = design)
fit = lmFit(datamatrices[["random_raw"]], design)
fit2 = contrasts.fit(fit, cont.matrix)
limma_p_rand_batchblocked = eBayes(fit2)$p.value[, 1]
```



```r
par(mfrow=c(1, 1))
thiscolors = c("red", "blue", "black")
hist(limma_p_rand_combat, breaks=100, freq=T, border=thiscolors[1], main="P-values, Random numbers for DP vs. N")
hist(limma_p_rand_nocombat, breaks=100, add=T, freq=T, border=thiscolors[2])
hist(limma_p_rand_batchblocked, breaks=100, add=T, freq=T, border=thiscolors[3])
legend("topright", legend=c("ComBat adjusted", "Limma adjusted", "No adjustment"), text.col=thiscolors)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 


The p-value plots from the ComBat adjusted random data shows a skewed distribution, although not as much as for the real data. The random data has no batch effect or other effects like the real data which might affect the ComBat adjustment.

Another sanity check is by permuting the DP (GA) and N(generic) labels and observe changes in the results. In order to retain the batch effect we will only swap labels within a batch. This should in theory give fewer significant probes.


```r
source("../../commonscripts/helperfunctions.r")

nshuffleddatasets = 10
datamatrices_permuted = list()
changedlabels = vector()
for (i in 1:nshuffleddatasets) {
    x = shufflesamplesinbatch(sampleannotation$code, sampleannotation$chip, 
        sampleannotation$covariate, c("N", "DP"))
    changedlabels[i] = (sum(sampleannotation[, "covariate"] != sampleannotation[x, 
        "covariate"]))  # see how many got change. Has 11 N's)  
    datamatrices_permuted[[i]] = as.matrix(ComBat(dat = datamatrices[["real_qnorm"]][, 
        x], batch = sampleannotation$chip, mod = combatmod, numCovs = NULL, 
        par.prior = TRUE, prior.plots = FALSE))
}
```

**10** permutations were created and batch adjusted with ComBat.

Then perform the limma test on the ComBat adjusted data as was done for the real data.

```r
permuted_pvalcounts = list()
overviewtab = data.frame(permutation = 1:nshuffleddatasets, changedlabels = changedlabels, 
    significantprobes = NA)
for (i in 1:nshuffleddatasets) {
    group = factor(sampleannotation$covariate)
    design = model.matrix(~0 + group)
    fit = lmFit(datamatrices_permuted[[i]], design)
    cont.matrix = makeContrasts(contrasts = "groupDP-groupN", levels = design)
    fit2 = contrasts.fit(fit, cont.matrix)
    this_p = eBayes(fit2)$p.value
    permuted_pvalcounts[[i]] = hist(this_p, plot = FALSE, breaks = 100)$counts
    overviewtab$significantprobes[i] = sum(p.adjust(this_p, method = "fdr") < 
        0.05)
}
print(overviewtab)
```

```
##    permutation changedlabels significantprobes
## 1            1             6               563
## 2            2            10               455
## 3            3             8               254
## 4            4             8               649
## 5            5             8               835
## 6            6             8               405
## 7            7             4               604
## 8            8            10              1130
## 9            9             8               238
## 10          10            12               662
```

```r
rm(overviewtab, group, design, fit, cont.matrix, fit2, this_p)
```

The table lists the number of labels that are changed compared to the real ones for each permutation. Remember that the total number of N (generic) samples is **11**. So even when about half of those are swapped with DP(GA), the ComBat adjusted data produces many significant genes. The number of significant genes using the same 0.05 FDR cut-off for the real labels was **1744**.

The P-value distribution can also be plotted.

```r
real_pvalcounts=hist(limma_p_alt, plot=FALSE, breaks=100)$counts
plot((1:20)/100, real_pvalcounts[1:20], col="red",
     main="P-values, permuted real data", 
      xlab="p-value", ylab="Frequency", type="l", lwd=2, 
      ylim=c(0, max(unlist( c(permuted_pvalcounts,real_pvalcounts)))))
for(i in 1:length(permuted_pvalcounts))
{
  lines((1:20)/100,permuted_pvalcounts[[i]][1:20], col="blue")
}
legend("topright", legend=c("Reproduced p-values real labels", "Permuted DP and N labels"), text.col=c("red", "blue"))
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27.svg) 

Shuffling the N and DP labels with in batches produces about the same p-value distribution as the real data. This is an indication that the significant genes found in the real data is not associated with the N and DP labels, but a rather an effect of the ComBat adjustments.



### References


Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127.

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
Platform: x86_64-w64-mingw32/x64 (64-bit)

locale:
[1] LC_COLLATE=English_United States.1252 
[2] LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] statmod_1.4.19        preprocessCore_1.24.0 limma_3.18.12        
 [4] sva_3.8.0             mgcv_1.7-28           nlme_3.1-113         
 [7] corpcor_1.6.6         GEOquery_2.28.0       Biobase_2.22.0       
[10] BiocGenerics_0.8.0    knitr_1.5            

loaded via a namespace (and not attached):
[1] evaluate_0.5.1  formatR_0.10    grid_3.0.2      lattice_0.20-24
[5] Matrix_1.1-2    RCurl_1.95-4.1  stringr_0.6.2   tools_3.0.2    
[9] XML_3.98-1.1   
```


generation ended 2014-05-24 00:14:56. Time spent 55 minutes .


