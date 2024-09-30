library(loopcolcox)
ica$OS_STATUS<-as.numeric(ica$OS_STATUS)

df<-coxbycol(ica$OS_MONTHS ,ica$OS_STATUS ,data)
head(df)

df%>%filter(significance=="YES" & prognosis == "unfavorable")->df2

## font.size : size of the font in the graph
coxvolcano(df,font.size=18)

df<-df[complete.cases(df),]
df%>%filter(coef.beta> -20)->df

library(patchwork)
## nb : number of covariates to put on the graph
p1<-plotbeta(df2,nb=16,title="",size=16)
p2<-plotnlphr(df2,nb=16,title="",size=16)
p1+p2