library(logitloop)
library(dplyr)

df<-logitloop(trans,outcome="group")
df

plotcoef(df,nb=20,title="")

write.table(df,file="logit.tsv",row.names=F,sep="\t")

df%>%filter(significance=="YES" & risk=="POSITIVE")->sig

##extract info for equation
id<-sig$predictors
beta<-sig$beta

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation



(ILMN_1719392_FH*4.469010534886)+(ILMN_1811367_MAT2B*4.25330754027876)+(ILMN_2410924_PLOD2*2.14866159055579)+
(ILMN_1684391_PLOD1*5.67387793796387)+(ILMN_1799139_PLOD2*2.75456081023585)+(ILMN_1790680_PDE6D*2.87447505556883)+
(ILMN_1755974_ALDOC*0.635624151457243)+(ILMN_1774281_NT5DC3*1.30199732330167)


library(multirocauc)
roc.list<-roclist(trans,id,outcome="group")
roc.list
rocplot(roc.list,line=1,title="metabolism ~ ICA - proliferative",police=14)

trans%>%mutate(meta.score=(ILMN_1719392_FH*4.469010534886)+(ILMN_1811367_MAT2B*4.25330754027876)+(ILMN_2410924_PLOD2*2.14866159055579)+
(ILMN_1684391_PLOD1*5.67387793796387)+(ILMN_1799139_PLOD2*2.75456081023585)+(ILMN_1790680_PDE6D*2.87447505556883)+
(ILMN_1755974_ALDOC*0.635624151457243)+(ILMN_1774281_NT5DC3*1.30199732330167))->trans

save(trans,file="metascore.rda")

library(ggplot2)
library(ggbeeswarm)

ggplot(trans,aes(group,meta.score))+geom_boxplot(outlier.shape=NA) + 
  scale_fill_brewer(palette="Dark2")+
  geom_point(aes(fill=factor(group),size=1),shape = 21, alpha = .8, position = position_dodge2(width = .5))+
  theme_classic(base_size=18) +
  theme(legend.position = "right")+xlab("ICA subclasses")+ggtitle("")