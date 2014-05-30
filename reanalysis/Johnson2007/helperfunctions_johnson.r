


loadjohnsondata = function()
{  
  datamatrix = as.matrix(read.table("data/dataExample2.txt", sep="\t", header=TRUE))
  sampleannotation = read.table("data/sampleInfoExample2.txt", 
                                sep="\t", header=TRUE, stringsAsFactors=FALSE)
  rownames(sampleannotation)=sampleannotation$ArrayName
  sampleannotation$Batch=factor(as.character(sampleannotation$Batch)) 
  datamatrix=datamatrix[,sampleannotation$Type!="WT"]
  sampleannotation=sampleannotation[sampleannotation$Type!="WT",]
  sampleannotation$Type=factor(sampleannotation$Type)
  
  log2data = datamatrix
  # flooring to 1
  log2data[log2data<1]=1
  # take out data with to much low/missing values.
  negativeprobesfilter =( rowSums(log2data>1) >= (0.9*ncol(log2data)) )
  log2data = log2data[negativeprobesfilter,]
  # quantilenormalize
  log2data=normalizeBetweenArrays(log2(log2data), method="quantile")
  
  return(list(sampleannotation=sampleannotation, rawdata=datamatrix, normdata=log2data))
}
