####GLMs to test effect of cage environment on OTU abundances at each age group
library("mvabund")
#Genus sequence abundance by elevation
generatable<-as.data.frame(read.csv(file.choose(), header = TRUE,row.names = 1))  #genera abundance by sample table
generameta <-read.csv(file.choose(), header=TRUE, row.names=1) #table of metadata for each sample (site and other environmental)
generaglm <- manyglm(abundance~group,data=generatable,family="neg binomial") #ran a glm with a negative binomial structure (there are tests you can do to see if your data fits the link function)
plot.manyglm(generaglm)
an.genera<-anova.manyglm(generaglm,test="wald",p.uni="adjusted")
an.genera
generadat <- mvabund(generatable)
plotMvaFactor(generadat, generameta$elevation)
