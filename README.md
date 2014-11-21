
README for the batch-adjust-warning-reports git repository.
----------------------


This git repository consists of more detailed case reports made in conjunction with, but
not part of, the article "Methods that remove batch effects while retaining group
differences may lead to exaggerated confidence in downstream analyses", V. Nygaard, E. A.
Rødland, E. Hovig, manuscript in preparation.

The target audience is the especially interested ones. The reports are not peer-reviewed and might be changed or updated after the publication data.

The format of the report is html, but the underlying r-code and text is available and
rerunnable. Run knit2html or the Rstudio knitr integration on the individual .rmd files to re-produce the plots and html.


### reanalysis/Johnson2007/

Reanalysis and alternative analysis of "Data set 2" from the original ComBat article:
Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray
expression data using Empirical Bayes methods. Biostatistics 8(1):118-127.

### reanalysis/Towfic2014/

Reanalysis and alternative analysis of parts of Towfic et al.(2014, January). Comparing
the biological impact of glatiramer acetate with the biological impact of a generic. PloS
one 9(1), e83757.

### reanalysis/Leek2012/

Full random number re-run of the tutorial from the sva package (v3.8.0), Leek JT,
Johnson WE, Parker HS, Jaffe AE, Storey JD.(2012) The sva package for removing batch
effects and other unwanted variation in high-throughput experiments. Bioinformatics. 2012
Mar 15;28(6):882-3.

### commonscripts/

- helperfunctions.r  ad hoc methods.

### Additional information

The scripts that produces the figures in "Methods that remove batch effects while retaining group
differences may lead to exaggerated confidence in downstream analyses" are found here:
https://github.com/ous-uio-bioinfo-core/batch-adjust-warning-figures.git
