


# Towfic et al.
# load data and probe annotation and some formatting
loadtowfic = function(downloaddata=TRUE)
{
  
  sampleannotation = read.table("data/sampleannotation.csv", sep="\t",
                                header=TRUE,  stringsAsFactors=FALSE)
  sampleannotation$code = make.names(sampleannotation$code)
  sampleannotation$chip = as.character(sampleannotation$chip)
  dimnames(sampleannotation)[[1]] = sampleannotation$code
  # take out 3 samples that are not assign to a geoaccession. Failed QC?
  sampleannotation = sampleannotation[!is.na(sampleannotation$geoaccession),] 
  
  
  if(downloaddata)
  {
    temp = tempfile()
    download.file(url="http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz",
                  destfile=temp, mode = "wb")
    rawdata = read.table(temp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    unlink(temp)
    
    temp = tempfile()
    download.file(url="http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file",
                  destfile=temp, mode = "wb")
    tardirtemp = tempfile()
    untar(temp, exdir = tardirtemp)
    rawannotation = read.table(paste(tardirtemp ,"/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt.gz", sep=""), 
                               sep="\t", header=TRUE, stringsAsFactors=FALSE, 
                               skip=8, comment.char="", quote="", fill=TRUE)
    
    unlink(temp)
    unlink(tardirtemp, recursive=TRUE)
  }else{
    
    # if download did not work change downloaddata to FALSE  
    # download data from the GEO deposit
    # unpack and place files in a folder named "not_in_github"
    rawdata = read.table("not_in_github/GSE40566_non_normalized.txt", 
                         sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # the probe annotation file found inside GSE40566_RAW.tar
    rawannotation = read.table("not_in_github/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt", 
                               sep="\t", header=TRUE, stringsAsFactors=FALSE, 
                               skip=8, comment.char="", quote="", fill=TRUE)
  }  
  
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
  
  return(list(sampleannotation=sampleannotation, data=datamatrix_raw, probes=probeannotation))
}