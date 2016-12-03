argv <- commandArgs(TRUE)
filedir<-argv[1]
files<-list.files(filedir)
require(minfi)
require(RPMM)
require(data.table)
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

bt<-matrix(rep(NA,485577*836),nrow=485577)
probe<-read.table("~/Data/BRCA/methy450k/tumor/jhu-usc.edu_BRCA.HumanMethylation450.9.lvl-3.TCGA-BH-A1F0-01A-11D-A138-05.txt",header=T,sep='\t')
probe<-probe[,1]
probe<-probe[-1]
rownames(bt)<-probe
require(data.table)
require(limma)
path_tumor<-"~/BRCAme450k/tumor"
path_normal<-"~/BRCAme450k/normal"
files_tumor<-list.files(path_tumor)
files_normal<-list.files(path_normal)
filespath_tumor<-paste(path_tumor,files_tumor,sep='/')
filespath_normal<-paste(path_normal,files_normal,sep='/')
filespath<-c(filespath_tumor,filespath_normal)
for(i in 1:length(filespath)){
  temp<-fread(filespath[i],sep='\t',skip=1)
  bt[temp[,V1],i]<-temp[,V4]
  rm(temp)
  print(i)
  gc()
}

num<-apply(bt,1,function(x) length(which(is.na(x)=="FALSE")))
bt<-bt[num==836,]

type<-c(rep("tumor",740),rep("normal",96))
design<-model.matrix(~factor(type))
colnames(design)<-c("Intercept","tumor")

fit<-lmFit(bt,design)
fit<-eBayes(fit)

design1<-design[,2]
fit1<-lmFit(bt,design1)
fit1<-eBayes(fit1)













