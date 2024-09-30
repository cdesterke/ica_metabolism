library(GEOquery)
# download data GEO GSE32225
gse <- getGEO("GSE32225", GSEMatrix = TRUE)

# if dataset has several data you can access like
set <- gse[[1]]

# extract expression
data<-exprs(set)
data<-log(data,2)
# extract annotation
features<-fData(set)

library(dplyr)
features%>%select(ID,Symbol)%>%distinct()->annot



# extract phenotype
pheno <- pData(set)
library(janitor)
pheno%>%clean_names()->pheno
library(skimr)
skim(pheno)


write.csv(pheno,file="pheno.csv")
pheno<-read.csv("pheno.csv",h=T,row.names=1)



# export data
all<-merge(annot,data,by="row.names")

write.csv(all,file="data.csv",row.names=F)


meta<-read.table("meta.tsv",h=T,sep="\t")

vector<-as.vector(meta$identifiers)

all%>%filter(Symbol%in%vector)->df

df$id<-paste(df$ID,df$Symbol,sep="_")

df%>%select(-c(1:3))%>%relocate(id)->df

row.names(df)<-df$id
df$id<-NULL

pheno$sub<-paste(pheno$tissue,pheno$group,sep="_")

all(row.names(pheno)==colnames(df))
library(transpipe)
pcatrans(df,pheno,group="sub")

bestheat(df,pheno,rownames=T,font=12)