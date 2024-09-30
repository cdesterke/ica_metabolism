library(loopcolcox)


df<-coxbycol(ica$OS_MONTHS ,ica$OS_STATUS ,data)
head(df)



## os score

os<-read.table(file="sigbadOS.tsv",h=T,sep="\t")




##extract info for equation
id<-df2$identifiers
beta<-df2$coef.beta

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation




library(dplyr)
## paste equation for quantitative variables 
data%>%mutate(metaos.score=(ACSM4*2.6030371391095)+(ALDOC*0.7228807445177)+(CASP9*1.70837919347801)+
(PLOD2*0.5676153037709)+(STEAP1B*0.413497016163516)+(PDE6D*2.18106351205589)+
(ERI1*2.28956138400499)+(ALPP*0.731410448276306)+(PNMT*0.648854679565766)+
(NT5M*0.573983311723712)+(NT5DC3*1.01227494713823)+(MAT2B*1.78532294973716)+
(GSTT2*0.367303156198813)+(PLOD1*0.91405270836529)+(FH*0.85166424842874)+
(AASS*0.508920750881935))->data

ica$metaos.score<-data$metaos.score
data$OS_MONTHS<-ica$OS_MONTHS
data$OS_STATUS<-ica$OS_STATUS
library(survminer)
res.cut <- surv_cutpoint(data, time = "OS_MONTHS", event = "OS_STATUS",variables = c("metaos.score"))
plot(res.cut, "metaos.score", palette = "Set1")
res.cat <- surv_categorize(res.cut)
library("survival")

fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~metaos.score, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE,palette="Set1",pval = TRUE,
risk.table.y.text=F,ggtheme=theme_classic2(base_size=16),conf.int.style="step")

ica$os.cat<-res.cat$metaos.score

ica$os.cat<-as.factor(ica$os.cat)
ica$os.cat<-relevel(ica$os.cat,ref="low")

save(ica,file="icametaosscore.rda")