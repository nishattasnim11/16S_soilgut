##GC data re-visualized with grouped barplot
##here I'm making my dataframe from scratch, but you can upload a csv and format it alternatively
SCFA=c("Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Acetic acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Butyric acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid",
       "Propionic acid"
       )

weight=c(0.10435637,
         0.19542415,
         0.22021168,
         0.05291482,
         0.05393569,
         0.08929114,
         0.09898163,
         0.06903061,
         0.06300686,
         0.14320670,
         0.12150779,
         0.10256354,
         
         0.029562248,
         0.060314598,
         0.093722865,
         0.052268921,
         0.069524964,
         0.005822619,
         0.084407127,
         0.229333741,
         0.016342242,
         0.191528252,
         0.021948484,
         0.131608980,
         
         0.032851537,
         0.115900230,
         0.134641393,
         0.028268984,
         0.037154284,
         0.061113052,
         0.045711007,
         0.028202721,
         0.035673404,
         0.080484140,
         0.067115396,
         0.058997240,
         
         0.018277685,
         0.029318160,
         0.054011280,
         0.021838374,
         0.029633223,
         0.003357772,
         0.064938254,
         0.121031426,
         0.002001930,
         0.106392770,
         0.007736106,
         0.041238869,
         
         0.016038638,
         0.023196390,
         0.024986103,
         0.009596128,
         0.011921697,
         0.011195844,
         0.014629573,
         0.005213550,
         0.006839255,
         0.020657433,
         0.016207873,
         0.014522346,
         
         0.008965786,
         0.013275355,
         0.016866904,
         0.010869425,
         0.009885849,
         0.004202399,
         0.018638153,
         0.033245206,
         0.003303628,
         0.028777121,
         0.002519504,
         0.014843349)

weeks=c(3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6)

df=data.frame(SCFA,weight)

##plotting grouped bar chart of mean %weight value with SD error bars 
df1  = transform(Forest, mean=rowMeans(df[cols]), sd=apply(df[cols],1, sd))

ggplot(df1, aes(x=as.factor(Gene), y=mean, fill=Species)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))

Forest = data.frame(SCFA, weight, weeks)
Forestplot<- ggplot(Forest, aes(factor(SCFA),weight)) + 
    geom_boxplot(aes(fill=factor(weeks)), show.legend=T)+ylab("Cecal SCFAs (% weight)")+
    xlab("")+
    guides(fill=guide_legend(title="Age(weeks)"))+
    geom_jitter()+
    theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
    

##doing a mann-whitney test on the null hypothesis that %weight cecal SCFA of 3 and 6 samples are identical populations 
#need to subset the dataset first
ForestAcetic<-subset(Forest, SCFA=="Acetic acid")
ForestButyric<-subset(Forest, SCFA=="Butyric acid")
ForestPropionic<-subset(Forest, SCFA=="Propionic acid")

#perform tests, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(weight ~ weeks, data=ForestAcetic) ##p-value 0.1432 so fail to reject null
wilcox.test(weight ~ weeks, data=ForestButyric) ##p-value 0.1135
wilcox.test(weight ~ weeks, data=ForestPropionic) ##p-value 0.5899


###Plot forest, urban and no soil for each SCFA
#need to do a bit of wrangling first
colnames(propionic3data)<-c("weight", "treatment")
colnames(propionic6data)<-c("weight", "treatment")
propionic3data$age<-3 ##add column "age" to dataframe and fill with value 3
propionic6data$age<-6 ##add column "age" to dataframe and fill with value 3
acetic3data$Age <- NULL ##delete column in dataframe by setting to NULL
Acetic <- rbind(acetic3data,acetic6data)
Butyric <- rbind(butyric3data,butyric6data)
Propionic <- rbind(propionic3data,propionic6data)

Aceticplot<- ggplot(Acetic, aes(factor(treatment),weight)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("Cecal SCFAs (% weight)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("Acetic Acid")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))

Butyricplot<- ggplot(Butyric, aes(factor(treatment),weight)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("Cecal SCFAs (% weight)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("Butyric Acid")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
Butyricplot

Propionicplot<- ggplot(Propionic, aes(factor(treatment),weight)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("Cecal SCFAs (% weight)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("Propionic Acid")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
Propionicplot

##subset
AceticUrban<-subset(Acetic, treatment=="urban")
AceticForest<-subset(Acetic, treatment=="forest")
AceticNS<-subset(Acetic, treatment=="no soil")

ButyricUrban<-subset(Butyric, treatment=="urban")
ButyricForest<-subset(Butyric, treatment=="forest")
ButyricNS<-subset(Butyric, treatment=="no soil")

PropionicUrban<-subset(Propionic, treatment=="urban")
PropionicForest<-subset(Propionic, treatment=="forest")
PropionicNS<-subset(Propionic, treatment=="no soil")

##test
wilcox.test(weight ~ age, data=AceticUrban) ##0.05158
wilcox.test(weight ~ age, data=AceticForest) ##0.1432 
wilcox.test(weight ~ age, data=AceticNS) ##0.7128

wilcox.test(weight ~ age, data=ButyricUrban) ##0.05156
wilcox.test(weight ~ age, data=ButyricForest) ##0.1135 
wilcox.test(weight ~ age, data=ButyricNS) ##0.7925

wilcox.test(weight ~ age, data=PropionicUrban) ##0.002914, reject null!!
wilcox.test(weight ~ age, data=PropionicForest) ##0.5899 
wilcox.test(weight ~ age, data=PropionicNS) ##1


##write out plot
pdffile <- "SCFA_AceticbyAge.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=10, height=5)
Aceticplot
dev.off()

pdffile <- "SCFA_ButyricbyAge.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=10, height=5)
Butyricplot
dev.off()

pdffile <- "SCFA_PropionicbyAge.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=10, height=5)
Propionicplot
dev.off()

pdffile <- "GC_Nov20.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
multiplot(Aceticplot,Propionicplot, Butyricplot,cols = 1)
dev.off()





