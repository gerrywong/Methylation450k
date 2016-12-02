argv <- commandArgs(TRUE)
filedir<-argv[1]
files<-list.files(filedir)
require(minfi)
require(RPMM)
source("BMIQ_1.3.R")
for(i in 1:length(files)){
  file<-files[i]
  sample<-strsplit(file,split="\\.")[[1]][6]
  temp<-readTCGA(paste(filedir,file,sep="/"))
  anno<-getAnnotation(temp)
  beta<-getBeta(temp)
  d<-data.frame(beta,type=as.numeric(as.factor(anno[,9])))
  d<-d[!is.na(d[,1]),]
  beta.v<-d[,1]
  design.v<-d[,2]
  result<-BMIQ(beta.v,design.v,plots=F,sampleID = sample)
  beta_BMIQ<-result[[1]]
  towrite<-cbind(d,beta_BMIQ)
  write.table(towrite,file=paste0(sample,"_BMIQ.txt"),quote=F,sep='\t')
}
