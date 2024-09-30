## multivariate
library(survival)
m<-coxph(formula=Surv(OS_MONTHS, OS_STATUS)~ISHAK_FIBROSIS_SCORE+AJCC_PATHOLOGIC_TUMOR_STAGE+metaos.score,data=ica)
m
summary(m)

library(broom)
library(broom.helpers)
library(GGally)


test<-cox.zph(m)
test
ggcoxzph(test,font.main = 8,ggtheme = theme_classic2(base_size=8))

ggcoef_model(m,exponentiate=T)+scale_color_brewer(palette="Dark2")+
	theme(text=element_text(size=14),legend.position="none")
library(rms)
ica%>%select(OS_MONTHS, OS_STATUS,metaos.score,ISHAK_FIBROSIS_SCORE,AJCC_PATHOLOGIC_TUMOR_STAGE)->small
ddist <- datadist(small)
oldoption <- options(datadist='ddist')


f<-cph(formula=Surv(OS_MONTHS, OS_STATUS)~metaos.score+ISHAK_FIBROSIS_SCORE+AJCC_PATHOLOGIC_TUMOR_STAGE,x=TRUE,y=TRUE,surv=TRUE,data=small)
surv<-Survival(f)

nomo <- nomogram(f,lp=TRUE,fun=function(x) surv(20,x),funlabel="20-Months Survival Prob")



plot(nomo)

set.seed(136879)
##m number of patient by groups, u=1 time 1 years
cal<-calibrate(f,B=500,method="boot",cmethod="KM",m=8,u=18)
plot(cal,xlab="Predicted probability at 18 months",ylab="Actual 18 months OS proportion" )

cal<-calibrate(f,B=500,method="boot",cmethod="KM",m=10,u=24)
plot(cal,xlab="Predicted probability at 2 year",ylab="Actual 24 months OS proportion" )



nomo <- nomogram(f,lp=TRUE,fun=list(function(x) surv(18,x)), funlabel=c("18 months Survival Prob"))

plot(nomo)
ggplot(Predict(f,metaos.score))+theme_minimal(base_size=16)+ggtitle("TCGA risk.score prediction")
library(regplot)
regplot(f, rank="sd",failtime = c(18))

res<-tidy(m,exponantial=T,conf.int = TRUE)

write.table(res,file="multivariatemetos.tsv",sep="\t",row.names=F)