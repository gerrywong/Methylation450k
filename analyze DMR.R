#####################################
#read BRCA tumor methylation450 file one by one
#because of out of memory by reading all files at once
#please setwd() in BRCA directory
#time 20160905
#by zxwang
#############################
library(data.table)
library(ggplot2)
library(minfi)
filepool<-list.files("BRCA_methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/")
filepath<-paste0("BRCA_methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/",filepool)
tem<-strsplit(filepool,split='\\.')
sample<-sapply(tem,"[",6)
ncol=length(filepool)
col_names<-sample
tumor_methy<-matrix(data=rep(0,485577*740),nrow=485577,ncol=ncol)
for(i in 1:ncol){
    temp<-fread(filepath[i],sep='\t',header=TRUE,skip=1)
    tumor_methy[,i]<-temp$Beta_value
    if(i==1){
        row_names<-temp$"Composite Element REF"
    }
    rm(temp)
    gc()
}
colnames(tumor_methy)<-col_names
rownames(tumor_methy)<-row_names
#tumor_methy_re<-tumor_methy[-which(is.na(tumor_methy[,1])==TRUE),]
#rm(tumor_methy)
gc()

filepool<-list.files("BRCA_methylation_normal/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/")
filepath<-paste0("BRCA_methylation_normal/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/",filepool)
tem<-strsplit(filepool,split='\\.')
sample<-sapply(tem,"[",6)
ncol=length(filepool)
col_names<-sample
normal_methy<-matrix(rep(0,485577*96),nrow=485577,ncol=96)
for(i in 1:ncol){
    temp<-fread(filepath[i],sep='\t',header=TRUE,skip=1)
    normal_methy[,i]<-temp$Beta_value
    if(i==1){
        row_names<-temp$"Composite Element REF"
    }
    rm(temp)
    gc()
}
rownames(normal_methy)<-row_names
colnames(normal_methy)<-col_names
#normal_methy_re<-normal_methy[-which(is.na(normal_methy[,1])==TRUE),]

normal_mean<-rowMeans(normal_methy)
normal_var<-apply(normal_methy,1,var)
normal_data<-cbind(normal_mean,normal_var)

tumor_mean<-rowMeans(tumor_methy)
tumor_var<-apply(tumor_methy,1,var)
##################
#paired z-test 90Vs90
sample_tumor<-substr(colnames(tumor_methy),1,12)
sample_normal<-substr(colnames(normal_methy),1,12)

paired<-intersect(sample_normal,sample_tumor)
choose_tumor<-which(sample_tumor %in% paired)
choose_normal<-which(sample_normal %in% paired)
#there is a problem: two tumor samples are from a same patient
#first method
choose_tumor1<-choose_tumor[-2]
choose_tumor2<-choose_tumor[-1]
#sort by name to make sure the columns with same order in two matrix
choose_normal_ordered<-choose_normal[order(sample_normal[choose_normal])]
choose_tumor1_ordered<-choose_tumor1[order(sample_tumor[choose_tumor1])]
choose_tumor2_ordered<-choose_tumor2[order(sample_tumor[choose_tumor2])]

#function
diff1<-normal_methy[,choose_normal_ordered]-tumor_methy[,choose_tumor1_ordered]
diff1_mean<-rowMeans(diff1)
diff1_sd<-apply(diff1,1,sd)
diff1_se<-diff1_sd/sqrt(ncol(diff1))
diff1_z<-diff1_mean/diff1_se
p1<-2*pnorm(abs(diff1_z),lower.tail=F)
diff1_zp<-cbind(diff1_z,p1)
diff1_zp_sort<-diff1_zp[order(diff1_zp[,2]),]
p1_adj<-p.adjust(diff1_zp_sort[,2],method="bonferroni")
diff1_zp<-cbind(diff1_zp_sort,p1_adj)
write.table(diff1_zp,file="BRCA_methy450_90Vs90_paired_ztest.txt",sep='\t',row.names=TRUE,quote=FALSE)
rm(diff1_zp_sort)
sig_probe1<-rownames(diff1_zp)[which(p1_adj<0.01)]

diff2<-normal_methy[,choose_normal_ordered]-tumor_methy[,choose_tumor2_ordered]
diff2_mean<-rowMeans(diff2)
diff2_sd<-apply(diff2,1,sd)
diff2_se<-diff2_sd/sqrt(ncol(diff2))
diff2_z<-diff2_mean/diff2_se
p2<-2*pnorm(abs(diff2_z),lower.tail=F)
diff2_zp<-cbind(diff2_z,p2)
diff2_zp_sort<-diff2_zp[order(diff2_zp[,2]),]
p2_adj<-p.adjust(diff2_zp_sort[,2],method="bonferroni")
diff2_zp<-cbind(diff2_zp_sort,p2_adj)
rm(diff2_zp_sort)

##96Vs650(649) z-test
#
indep_tumor_methy<-tumor_methy[,-choose_tumor]
indep_tumor_mean<-rowMeans(indep_tumor_methy)
indep_tumor_var<-apply(indep_tumor_methy,1,var)

n1<-ncol(normal_methy)
n2<-ncol(indep_tumor_methy)
z_test<-(normal_mean-indep_tumor_mean)/sqrt(normal_var/n1+indep_tumor_var/n2)
p<-2*pnorm(abs(z_test),lower.tail=F)
p_adj<-p.adjust(p,method="bonferroni")
result<-cbind(z_test,p,p_adj)
result<-result[order(result[,2]),]
write.table(result,file="BRCA_methy450_96Vs649_test.txt",sep='\t',row.names=TRUE,quote=FALSE)

#get the signature probes responding genes and chromsome
temp<-fread("jhu_**.txt",sep='\t',header=TRUE,skip=1)
colnames(temp)<-c("Composite_Element_REF","Beta_value","Gene_Symbol","Chromosome","Genomic_Coordinate")
sig_probe2<-rownames(diff2_zp)[which(p2_adj<0.01)]
#loc<-match(sig_probe2,temp[,get(names(temp)[1])])
setkey(temp,"Composite_Element_REF")
sig_gene<-temp[sig_probe2,.(Composite_Element_REF,Gene_Symbol,Chromosome,Genomic_Coordinate)]#注意这里取列的方式
sig_gene<-sig_gene[order(sig_gene$Chromosome,sig_gene$Genomic_Coordinate),]

sig_table<-table(sig_gene$Gene_Symbol)
choose_gene<-names(sig_table)[which(sig_table>4)]
choose_probe<-sig_gene[Gene_Symbol %in% choose_gene]$Composite_Element_REF
cat(choose_probe,file="Choose_probe_of_90_paired_test.txt",sep='\t')

pdf(file="heatmap_of_paired_test.pdf",width=9,height=7)
heatmap.2(normal_methy[choose_probe,])
dev.off()

##validate that the purity can affect the detection of dmp
purity<-read.table("BRCA-purity.txt",header=T,sep='\t')
purity2<-purity[-which(!purity[,1]%in%colnames(tumor_methy)),]
purity2<-purity2[order(purity2[,2]),]
topten<-tumor_methy[,as.character(purity2[667:740,1])]
bottomten<-tumor_methy[,as.character(purity2[1:74,1])]
normalmethy2<-normal_methy[,-c(39,40,41,42,44,45,46,48)] ## 86

##直接使用默认计算均值的方法，有NA值的就不计算
z_test<-function(normal,tumor){
  normal_mean<-rowMeans(normal)
  normal_var<-apply(normal,1,var)
  tumor_mean<-rowMeans(tumor)
  tumor_var<-apply(tumor,1,var)
  ztest<-(normal_mean-tumor_mean)/sqrt(normal_var/ncol(normal)+tumor_var/ncol(tumor))
  p<-2*pnorm(abs(ztest),lower.tail=F)
  p_adj<-p.adjust(p,method="bonferroni")
  res<-cbind(ztest,p_adj)
  return(res)
}

###topten Vs normalmethy2
toptenztest<-z_test(normalmethy2,topten)
###bottomten Vs normalmethy2
bottomtenztest<-z_test(normalmethy2,bottomten)

toptenztest<-toptenztest[order(toptenztest[,2]),]
bottomtenztest<-bottomtenztest[order(bottomtenztest[,2]),]

toptendmp<-row.names(toptenztest[1:10000,])
bottomtendmp<-row.names(bottomtenztest[1:10000,])

library(minfi)
anno<-mapToGenome(file)
color<-rep("top",10000)
zvalue<-toptenztest[1:10000,1]
pvalue<-toptenztest[1:10000,2]
toptendmptoPlot<-data.frame(dmp=toptendmp,anno[toptendmp,c(1,2)],zvalue,pvalue,color)
colnames(toptendmptoPlot)<-c("dmp","chr","pos","zvalue","adjp","color")
color<-rep("bottom",10000)
zvalue<-bottomtenztest[1:10000,1]
pvalue<-bottomtenztest[1:10000,2]
bottomtentoPlot<-data.frame(dmp=bottomtendmp,anno[bottomtendmp,c(1,2)],zvalue,pvalue,color)
colnames(bottomtentoPlot)<-c("dmp","chr","pos","zvalue","adjp","color")
plotdata<-rbind(toptendmptoPlot,bottomtentoPlot,make.row.names=F)
colnames(plotdata)<-c("dmp","chr","pos","zvalue","adjp","color")
plotdata$chr<-as.factor(plotdata$chr)
library(ggplot2)
pdf(file="ztest by chromsome.pdf",width=12,height=12)
p<-ggplot(plotdata,aes(x=pos,y=zvalue,color=color))+geom_point(size=0.9)+facet_wrap(~chr,scales="free_x")
print(p)
dev.off()

common<-intersect(toptendmp,bottomtendmp)
subplotdata<-plotdata[which(plotdata$dmp%in%common),]
ggplot(subplotdata,aes(color,zvalue,fill=color))+geom_boxplot()+ggtitle("common1615 dmp zvalue")

##Mr. Teng  this plot can illstrate that purity of samples can affect the detection of differential methylation positions.
idx=!is.na(bottomtenztest[,1])&!is.na(toptenztest[,1])
plot(density(bottomtenztest[dix,1],na.rm=T),xlim=c(-40,40))
lines(density(toptenztest[idx,1],na.rm=T),col='red')
abline(v=c(-13,10),lty=2)

##plot the barplot of bottom&top samples
library(data.table)
toptentoplot<-melt(topten)
ggplot(toptentoplot,aes(value,y=..density..))+geom_histogram(binwidth = 0.005)
  +geom_density(kernel="gaussian")+facet_wrap(~Var2,scales="free_y")

bottomtentoplot<-melt(bottomten)
pdf("bottomten_samples_density.pdf",width=24,height=24)
p<-ggplot(bottomtentoplot,aes(value,y=..density..))+geom_histogram(binwidth = 0.005)
  +geom_density(kernel="gaussian")+facet_wrap(~Var2,scales="free_y")
print(p)
dev.off()

diffofztest<-toptenztestorderbyrowname[,1]-bottomtenztestorderbyrowname[,1]
cat(names(diffofztest[which(abs(diffofztest)>20),]),file="rownames_of_diff_more__than_20.txt")
probemorethantwt<-scan("/media/gerry/data2/Data/rownames_of_diff_more__than_20.txt",what="character")
bottomtentoplot<-bottomten[probemorethantwt,]
toptentoplot<-topten[probemorethantwt,]
bottomtentoplot<-data.frame(probe=rownames(bottomtentoplot),bottomtentoplot,rep("bottom",2138))
toptentoplot<-data.frame(probe=rownames(toptentoplot),toptentoplot,rep("top",2138))
colnames(bottomtentoplot)<-c("probe",1:74,"color")
colnames(toptentoplot)<-c("probe",1:74,"color")
mydata<-rbind(bottomtentoplot,toptentoplot)
mydata2<-melt(mydata,id.vars=c("probe","color"))
mydata2<-mydata2[order(mydata2[,1]),]
ggplot(mydata2[1:1480,],aes(probe,value,fill=color))+geom_boxplot()

##plot all normal samples' beta value of I&II probe
for(i in 1:ncol(normalmethy2)){
  toplot<-data.frame(value=normalmethy2[rownames(anno),i],type=anno[,9])
  colnames(toplot)<-c("value","type")
  ggplot(toplot,aes(value,color=type))+geom_density()+ggtitle(colnames(normalmethy2)[i])
  ggsave(file=paste0(colnames(normalmethy2)[i],".png"))
}
















