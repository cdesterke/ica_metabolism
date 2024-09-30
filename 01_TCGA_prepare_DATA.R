load("all.rda")
library(dplyr)

table(all$CANCER_TYPE_DETAILED)

all%>%filter(CANCER_TYPE_DETAILED=="Intrahepatic Cholangiocarcinoma")->ica



meta<-read.csv("metabolism.csv",h=T)

## biomart 112  gene annotation for metabolism database
bm<-read.table("bm112.txt", h=T,sep="\t",na.strings = "NA")


head(bm,n=100)

bm[bm==""] <- NA

bm<-bm[complete.cases(bm),]

dim(bm)

head(bm)

bm112<-bm

bm112%>%inner_join(meta,c("Mouse.gene.name"="mm.gene"))->metabolism

save(bm112,file="bm112converter.rda")

metabolism%>%select(Gene.name)%>%distinct()->sel
vector<-sel$Gene.name
inter<-intersect(vector,colnames(ica))

ica%>%select(all_of(inter))->data
data<-log(data+1,2)