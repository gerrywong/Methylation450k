#####################################
#read BRCA tumor methylation450 file one by one
#because of out of memory by reading all files at once
#please setwd() in BRCA directory
#time 20160905
#by zxwang
#############################
library("data.table")
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












