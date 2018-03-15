##Local cytokines with qPCR

reg3<-c(0.54150, 0.13798, 0.13732, 0.18839, 0.26390, 0.15347,0.24379,0.54083,0.20185,0.56382,1.34881,1.62006,0.14704,0.17123,0.26090)
reg3groups=rep(c("forest","no soil","urban"),times=c(6,3,6))
reg3data<-data.frame(reg3,reg3groups)
r<- ggplot(reg3data, aes(factor(reg3groups),reg3))
##increased font size of axis text and text to 16
r<-r+geom_boxplot(aes(fill=factor(reg3groups)),show.legend = F)+ylab("Relative expression")+xlab("")+ggtitle("Reg3-gamma")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
r
kruskal.test(reg3~reg3groups,reg3data)
reg3.mod1=lm(reg3~reg3groups, data=reg3data)
summary(reg3.mod1)
anova(reg3.mod1)

IFN<-c(0.31298,
       0.99165,
       0.90415,
       0.64201,
       1.43114,
       0.25169,
       0.42736,
       0.16776,0.76906,0.35188,0.31668,0.32089,0.24593,2.23666,1.62006,0.96354,1.08151,0.87700)
IFNgroups=rep(c("forest","no soil","urban"),times=c(8,3,7))
IFNdata<-data.frame(IFN,IFNgroups)
r2<- ggplot(IFNdata, aes(factor(IFNgroups),IFN))
r2<-r2+geom_boxplot(aes(fill=factor(IFNgroups)),show.legend = F)+ylab("Relative expression")+xlab("")+ggtitle("IFN-gamma")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
t.test(IFN~IFNgroups,IFNdata)


IL10<-c(0.22400, 1.27215, 1.76073,0.27303,0.55445,0.35671,0.30819,0.30935,0.93568,0.14409,0.89064,0.45427,0.97882,1.79600,1.62006,0.82850,1.27443,0.41095)
IL10data<-data.frame(IL10,IFNgroups)
r3<- ggplot(IL10data, aes(factor(IFNgroups),IL10))
r3<-r3+geom_boxplot(aes(fill=factor(IFNgroups)),show.legend = F)+ylab("Relative expression")+xlab("")+ggtitle("IL-10")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
kruskal.test(IL10~IFNgroups,IL10data)

TGFB<-c(1.65465,
        1.10444,
        1.54221,
        1.01759,
        0.84700,
        1.05665,
        1.12167,
        0.81429,1.29367,1.21788,1.67126,1.68006,1.17295,1.01627,0.99899,1.07082,1.10099,1.20827)
TGFBdata<-data.frame(TGFB,IFNgroups)
r4<- ggplot(TGFBdata, aes(factor(IFNgroups),TGFB))
r4<-r4+geom_boxplot(aes(fill=factor(IFNgroups)),show.legend = F)+ylab("Relative expression")+xlab("")+ggtitle("TGF-beta")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
kruskal.test(TGFB~IFNgroups,TGFBdata)

TNFa<-c(0.77857,
        1.86336,
        1.08042,
        0.89048,
        0.98470,
        1.11082,
        1.13575,
        0.43807,0.92877,0.87587,1.04167,1.55167,1.31282,0.81569,1.62006,0.40079,0.51432,1.36956)
TNFadata<-data.frame(TNFa,IFNgroups)
r5<- ggplot(TNFadata, aes(factor(IFNgroups),TNFa))
r5<-r5+geom_boxplot(aes(fill=factor(IFNgroups)),show.legend = F)+ylab("Relative expression")+xlab("")+ggtitle("TNF-alpha")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
kruskal.test(TNFa~IFNgroups,TNFadata)

IL1B<-c(0.70153,
        2.18415,
        1.09643,
        1.33861,
        0.51412,
        0.93362,
        0.64949,
        0.21873,0.71603,0.32206,0.39615,0.40373,0.74662,0.62385,0.50663,0.53698,0.37277,0.26721)
IL1Bdata<-data.frame(IL1B,IFNgroups)
r6<- ggplot(IL1Bdata, aes(factor(IFNgroups),IL1B))
r6<-r6+geom_boxplot(aes(fill=factor(IFNgroups)),show.legend = F)+ylab("Relative expression")+xlab("")+ggtitle("IL1B")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
kruskal.test(IL1B~IFNgroups,IL1Bdata)
boxplot(IL1B~IFNgroups,IL1Bdata)
multiplot(r,r2,r3,r4,r5,r6,cols = 3)

##Systemic cytokines
Eotaxin<-c(269.02,
           196.13,
           274.62,
           444.69,
           246.41,
           309.20,
           267.11,
           258.50,
           490.89,
           300.46,
           242.42,
           521.22,
           238.99,
           386.25,
           178.39,
           240.87,
           457.02,
           108.93,
           231.17,
           147.46,
           367.07,
           235.71,
           276.12,
           122.36)
Eotaxingroups=rep(c("no soil","urban","forest"),times=c(8,8,8))
Eotaxindata<-data.frame(Eotaxin,Eotaxingroups)
s<- ggplot(Eotaxindata, aes(factor(Eotaxingroups),Eotaxin))
s<-s+geom_boxplot(aes(fill=factor(Eotaxingroups)),show.legend = F)+ylab("Eotaxin (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s
kruskal.test(Eotaxin~Eotaxingroups,Eotaxindata)

Eotaxin6<-c(650.51,
            325.56,
            572.03,
            414.67,
            602.82,
            499.73,
            658.69,
            633.59,
            744.56,
            438.61,
            644.88,
            270.36,
            571.34,
            492.15,
            464.39,
            394.60,
            593.75,
            223.48,
            365.11)
Eotaxin6groups=rep(c("urban","forest","no soil"),times=c(8,8,3))
Eotaxin6data<-data.frame(Eotaxin6,Eotaxin6groups)
s61<- ggplot(Eotaxin6data, aes(factor(Eotaxin6groups),Eotaxin6))
s61<-s61+geom_boxplot(aes(fill=factor(Eotaxin6groups)),show.legend = F)+ylab("Eotaxin (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s61
kruskal.test(Eotaxin6~Eotaxin6groups,Eotaxin6data)

##wrangle
colnames(Eotaxindata)<-c("cytokine", "treatment")
colnames(Eotaxin6data)<-c("cytokine", "treatment")
Eotaxindata$age<-3
Eotaxin6data$age<-6
EotaxinTotal <- rbind(Eotaxindata,Eotaxin6data)

##plot 
EotaxinTotalplot<- ggplot(EotaxinTotal, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("Eotaxin")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
EotaxinTotalplot

##subset
EF<-subset(EotaxinTotal, treatment=="forest")
EU<-subset(EotaxinTotal, treatment=="urban")
EN<-subset(EotaxinTotal, treatment=="no soil")

##wilcoxon
wilcox.test(cytokine ~ age, data=EF) ##0.003 reject null!
wilcox.test(cytokine ~ age, data=EU) ##0.005 reject null!
wilcox.test(cytokine ~ age, data=EN) ##0.497

GCSF<-c(373.27,
        299.77,
        259.46,
        354.73,
        1117.08,
        314.40,
        345.25,
        564.68,
        643.95,
        567.38,
        576.56,
        482.87,
        453.75,
        412.04,
        343.88,
        485.16,
        522.43,
        473.64,
        512.37,
        414.51,
        742.88,
        296.04,
        549.99,
        379.77)
GCSFgroups=rep(c("no soil","urban","forest"),times=c(8,8,8))
GCSFdata<-data.frame(GCSF,GCSFgroups)
s2<- ggplot(GCSFdata, aes(factor(GCSFgroups),GCSF))
s2<-s2+geom_boxplot(aes(fill=factor(GCSFgroups)),show.legend = F)+ylab("G-CSF (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s2
kruskal.test(GCSF~GCSFgroups,GCSFdata)

GCSF6<-c(814.98,
         740.80,
         406.44,
         437.06,
         237.75,
         257.83,
         426.16,
         1394.56,
         1337.93,
         860.03,
         2264.58,
         695.61,
         509.00,
         781.42,
         574.94,
         520.20,
         653.39,
         153.82,
         464.33)
GCSF6groups=rep(c("urban","forest","no soil"),times=c(8,8,3))
GCSF6data<-data.frame(GCSF6,GCSF6groups)
s62<- ggplot(GCSF6data, aes(factor(GCSF6groups),GCSF6))
s62<-s61+geom_boxplot(aes(fill=factor(GCSF6groups)),show.legend = F)+ylab("G-CSF(pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s62
kruskal.test(GCSF6~GCSF6groups,GCSF6data)


##wrangle
colnames(GCSFdata)<-c("cytokine", "treatment")
colnames(GCSF6data)<-c("cytokine", "treatment")
GCSFdata$age<-3
GCSF6data$age<-6
GCSFTotal <- rbind(GCSFdata,GCSF6data)

##plot 
GCSFTotalplot<- ggplot(GCSFTotal, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("G-CSF")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
GCSFTotalplot

##subset
GCSFF<-subset(GCSFTotal, treatment=="forest")
GCSFU<-subset(GCSFTotal, treatment=="urban")
GCSFN<-subset(GCSFTotal, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=GCSFF) ## 0.015 reject null! 
wilcox.test(cytokine ~ age, data=GCSFU) ## 0.80
wilcox.test(cytokine ~ age, data=GCSFN) ## 0.92

IL1a<-c(235.44,
        210.13,
        284.79,
        271.36,
        342.11,
        316.99,
        284.79,
        301.45,
        246.54,
        453.07,
        195.05,
        331.63,
        81.12,
        331.63,
        412.90,
        455.63,
        153.05,
        153.05,
        256.90,
        177.06,
        266.67,
        261.85,
        324.41,
        195.05)
IL1agroups=rep(c("no soil","urban","forest"),times=c(8,8,8))
IL1adata<-data.frame(IL1a,IL1agroups)
s4<- ggplot(IL1adata, aes(factor(IL1agroups),IL1a))
s4<-s4+geom_boxplot(aes(fill=factor(IL1agroups)),show.legend = F)+ylab("IL-1a (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s4
kruskal.test(IL1a~IL1agroups,IL1adata)

IL1a6<-c(177.06,
         153.05,
         39.76,
         293.28,
         275.94,
         251.80,
         458.18,
         485.44,
         301.45,
         338.66,
         383.84,
         195.05,
         386.85,
         153.05,
         669.59,
         223.41,
         223.41,
         246.54,
         153.05)
IL1a6groups=rep(c("urban","forest","no soil"),times=c(8,8,3))
IL1a6data<-data.frame(IL1a6,IL1a6groups)
s64<- ggplot(IL1a6data, aes(factor(IL1a6groups),IL1a6))
s64<-s64+geom_boxplot(aes(fill=factor(IL1a6groups)),show.legend = F)+ylab("IL-1a (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s64
kruskal.test(IL1a6~IL1a6groups,IL1a6data)

##wrangle
colnames(IL1adata)<-c("cytokine", "treatment")
colnames(IL1a6data)<-c("cytokine", "treatment")
IL1adata$age<-3
IL1a6data$age<-6
IL1aTotal <- rbind(IL1adata,IL1a6data)

##plot 
IL1aTotalplot<- ggplot(IL1aTotal, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("IL-1a")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
IL1aTotalplot

##subset
IL1aF<-subset(IL1aTotal, treatment=="forest")
IL1aU<-subset(IL1aTotal, treatment=="urban")
IL1aN<-subset(IL1aTotal, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=IL1aF) ## 0.11 
wilcox.test(cytokine ~ age, data=IL1aU) ## 0.63
wilcox.test(cytokine ~ age, data=IL1aN) ## 0.08


IL1bS<-c(8.39,
        6.20,
        8.89,
        6.20,
        2.61,
        5.42,
        9.42,
        7.46,
        8.39,
        4.69,
        4.69,
        3.98,
        7.03,
        5.05,
        5.42,
        2.61,
        2.28,
        6.20,
        3.98,
        7.03,
        4.69,
        5.42,
        7.92,
        3.29)
IL1bSgroups=rep(c("no soil","urban","forest"),times=c(8,8,8))
IL1bSdata<-data.frame(IL1bS,IL1bSgroups)
s3<- ggplot(IL1bSdata, aes(factor(IL1bSgroups),IL1bS))
s3<-s3+geom_boxplot(aes(fill=factor(IL1bSgroups)),show.legend = F)+ylab("IL-1B (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s3
kruskal.test(IL1bS~IL1bSgroups,IL1bSdata)

IL1b6<-c(8.89,
         5.42,
         7.03,
         3.63,
         8.89,
         4.69,
         3.98,
         5.42,
         4.33,
         7.46,
         7.03,
         5.42,
         8.89,
         3.98,
         5.42,
         1.94,
         16.56)
IL1b6groups=rep(c("urban","forest","no soil"),times=c(6,8,3))
IL1b6data<-data.frame(IL1b6,IL1b6groups)
s63<- ggplot(IL1b6data, aes(factor(IL1b6groups),IL1b6))
s63<-s63+geom_boxplot(aes(fill=factor(IL1b6groups)),show.legend = F)+ylab("IL-1B (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s63
kruskal.test(IL1b6~IL1b6groups,IL1b6data)

##wrangle
colnames(IL1bSdata)<-c("cytokine", "treatment")
colnames(IL1b6data)<-c("cytokine", "treatment")
IL1bSdata$age<-3
IL1b6data$age<-6
IL1bTotal <- rbind(IL1bSdata,IL1b6data)

##plot 
IL1bTotalplot<- ggplot(IL1bTotal, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("IL-1B")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
IL1bTotalplot

##subset
IL1bF<-subset(IL1bTotal, treatment=="forest")
IL1bU<-subset(IL1bTotal, treatment=="urban")
IL1bN<-subset(IL1bTotal, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=IL1bF) ## 0.5
wilcox.test(cytokine ~ age, data=IL1bU) ## 0.3
wilcox.test(cytokine ~ age, data=IL1bN) ## 0.7

IL10S<-c(2.16,
         5.55,
         6.23,
         6.23,
         4.28,
         5.55,
         3.66,
         5.55,
         2.16,
         0.69,
         2.46,
         4.28,
         1.87,
         2.16,
         6.23,
         6.93,
         1.87,
         3.06,
         2.46,
         3.66,
         1.28,
         3.66,
         3.66)
IL10Sgroups=rep(c("no soil","urban","forest"),times=c(7,8,8))
IL10Sdata<-data.frame(IL10S,IL10Sgroups)
s5<- ggplot(IL10Sdata, aes(factor(IL10Sgroups),IL10S))
s5<-s5+geom_boxplot(aes(fill=factor(IL10Sgroups)),show.legend = F)+ylab("IL-10 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s5
kruskal.test(IL10S~IL10Sgroups,IL10Sdata)

IL106<-c(4.28,
         6.23,
         4.28,
         3.06,
         4.91,
         1.87,
         7.67,
         1.87,
         7.67,
         1.87,
         4.91,
         2.46,
         3.36,
         4.28,
         4.28,
         4.91)
IL106groups=rep(c("urban","forest","no soil"),times=c(6,7,3))
IL106data<-data.frame(IL106,IL106groups)
s65<- ggplot(IL106data, aes(factor(IL106groups),IL106))
s65<-s65+geom_boxplot(aes(fill=factor(IL106groups)),show.legend = F)+ylab("IL-10 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s65
kruskal.test(IL106~IL106groups,IL106data)

##wrangle
colnames(IL10Sdata)<-c("cytokine", "treatment")
colnames(IL106data)<-c("cytokine", "treatment")
IL10Sdata$age<-3
IL106data$age<-6
IL10Total <- rbind(IL10Sdata,IL106data)

##plot 
IL10Totalplot<- ggplot(IL10Total, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("IL-10")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
IL10Totalplot

##subset
IL10F<-subset(IL10Total, treatment=="forest")
IL10U<-subset(IL10Total, treatment=="urban")
IL10N<-subset(IL10Total, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=IL10F) ## 0.6
wilcox.test(cytokine ~ age, data=IL10U) ## 0.3
wilcox.test(cytokine ~ age, data=IL10N) ## 0.6

MCP1<-c(21.51,
       16.42,
       80.77,
       46.08,
       11.40,
       26.69,
       13.90,
       21.51,
       26.69,
       43.17,
       32.00,
       32.00,
       43.17,
       34.72,
       26.69,
       32.00,
       34.72,
       26.69,
       1.40,
       32.00,
       67.66,
       26.69,
       29.32,
       11.40)
MCP1groups=rep(c("no soil","urban","forest"),times=c(8,8,8))
MCP1data<-data.frame(MCP1,MCP1groups)
s6<- ggplot(MCP1data, aes(factor(MCP1groups),MCP1))
s6<-s6+geom_boxplot(aes(fill=factor(MCP1groups)),show.legend = F)+ylab("MCP-1 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s6
kruskal.test(MCP1~MCP1groups,MCP1data)

MCP16<-c(87.44,
         43.17,
         37.49,
         43.17,
         32.00,
         16.42,
         147.78,
         67.66,
         49.04,
         16.42,
         49.04,
         34.72,
         49.04,
         18.95,
         16.42,
         32.00,
         16.42,
         64.46)
MCP16groups=rep(c("urban","forest","no soil"),times=c(8,7,3))
MCP16data<-data.frame(MCP16,MCP16groups)
s66<- ggplot(MCP16data, aes(factor(MCP16groups),MCP16))
s66<-s66+geom_boxplot(aes(fill=factor(MCP16groups)),show.legend = F)+ylab("MCP-1 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s66
kruskal.test(MCP16~MCP16groups,MCP16data)

##wrangle
colnames(MCP1data)<-c("cytokine", "treatment")
colnames(MCP16data)<-c("cytokine", "treatment")
MCP1data$age<-3
MCP16data$age<-6
MCP1Total <- rbind(MCP1data,MCP16data)

##plot 
MCP1Totalplot<- ggplot(MCP1Total, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("MCP-1")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
MCP1Totalplot

##subset
MCP1F<-subset(MCP1Total, treatment=="forest")
MCP1U<-subset(MCP1Total, treatment=="urban")
MCP1N<-subset(MCP1Total, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=MCP1F) ## 0.5
wilcox.test(cytokine ~ age, data=MCP1U) ## 0.1
wilcox.test(cytokine ~ age, data=MCP1N) ## 0.5

Tnfa<-c(15.22,
        18.31,
        13.35,
        13.35,
        17.45,
        6.91,
        10.33,
        2.59,
        16.13,
        6.91,
        12.37,
        13.35,
        12.37,
        17.01,
        12.37,
        13.35,
        14.76,
        15.68,
        13.83,
        8.11,
        11.87,
        11.87,
        11.37,
        13.35)
Tnfagroups=rep(c("no soil","urban","forest"),times=c(8,8,8))
Tnfadata<-data.frame(Tnfa,Tnfagroups)
s7<- ggplot(Tnfadata, aes(factor(Tnfagroups),Tnfa))
s7<-s7+geom_boxplot(aes(fill=factor(Tnfagroups)),show.legend = F)+ylab("TNF-alpha (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s7
kruskal.test(Tnfa~Tnfagroups,Tnfadata)

Tnfa6<-c(12.37,
         10.33,
         15.22,
         11.37,
         8.11,
         5.63,
         17.01,
         12.86,
         15.22,
         13.35,
         10.85,
         13.35,
         6.91,
         6.91,
         2.59,
         6.91,
         4.94,
         8.11)
Tnfa6groups=rep(c("urban","forest","no soil"),times=c(7,8,3))
Tnfa6data<-data.frame(Tnfa6,Tnfa6groups)
s67<- ggplot(Tnfa6data, aes(factor(Tnfa6groups),Tnfa6))
s67<-s67+geom_boxplot(aes(fill=factor(Tnfa6groups)),show.legend = F)+ylab("TNF-alpha (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s67
kruskal.test(Tnfa6~Tnfa6groups,Tnfa6data)

##wrangle
colnames(Tnfadata)<-c("cytokine", "treatment")
colnames(Tnfa6data)<-c("cytokine", "treatment")
Tnfadata$age<-3
Tnfa6data$age<-6
TnfaTotal <- rbind(Tnfadata,Tnfa6data)

##plot 
TnfaTotalplot<- ggplot(TnfaTotal, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("TNF-alpha")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
TnfaTotalplot

##subset
TnfaF<-subset(TnfaTotal, treatment=="forest")
TnfaU<-subset(TnfaTotal, treatment=="urban")
TnfaN<-subset(TnfaTotal, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=TnfaF) ## 0.26
wilcox.test(cytokine ~ age, data=TnfaU) ## 0.32
wilcox.test(cytokine ~ age, data=TnfaN) ## 0.15

IL13<-c(16.76,
        33.95,
        24.61,
        28.59,
        31.26,
        21.98,
        29.92,
        33.28,
        42.14,
        29.92,
        28.59,
        21.98,
        35.31,
        20.67,
        33.95,
        33.95,
        33.95,
        32.61,
        16.76,
        23.29,
        28.59,
        35.31,
        25.27,
        32.61)
IL13groups=rep(c("no soil","urban","forest"),times=c(8,8,8))
IL13data<-data.frame(IL13,IL13groups)
s8<- ggplot(IL13data, aes(factor(IL13groups),IL13))
s8<-s8+geom_boxplot(aes(fill=factor(IL13groups)),show.legend = F)+ylab("IL-13 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s8
kruskal.test(IL13~IL13groups,IL13data)

IL136<-c(17.41,
         15.46,
         31.26,
         27.92,
         28.59,
         12.87,
         21.98,
         21.32,
         29.92,
         21.98,
         18.06,
         33.95,
         14.16,
         30.59,
         22.63,
         15.46,
         24.61,
         35.31,
         31.93)
IL136groups=rep(c("urban","forest","no soil"),times=c(8,8,3))
IL136data<-data.frame(IL136,IL136groups)
s68<- ggplot(IL136data, aes(factor(IL136groups),IL136))
s68<-s68+geom_boxplot(aes(fill=factor(IL136groups)),show.legend = F)+ylab("IL-13 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s68
kruskal.test(IL136~IL136groups,IL136data)

##wrangle
colnames(IL13data)<-c("cytokine", "treatment")
colnames(IL136data)<-c("cytokine", "treatment")
IL13data$age<-3
IL136data$age<-6
IL13Total <- rbind(IL13data,IL136data)

##plot 
IL13Totalplot<- ggplot(IL13Total, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("IL-13")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
IL13Totalplot

##subset
IL13F<-subset(IL13Total, treatment=="forest")
IL13U<-subset(IL13Total, treatment=="urban")
IL13N<-subset(IL13Total, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=IL13F) ## 0.14
wilcox.test(cytokine ~ age, data=IL13U) ## 0.03 reject null!
wilcox.test(cytokine ~ age, data=IL13N) ## 0.41

IL6<-c(2.63,
       2.70,
       3.73,
       3.89,
       8.45,
       4.05,
       4.21,
       3.42,
       6.12,
       4.69,
       3.42,
       4.77,
       2.31,
       3.89,
       2.63,
       4.37,
       3.26,
       4.21,
       2.63,
       6.76,
       5.00,
       2.63,
       5.08,
       14.97)
IL6groups=rep(c("no soil","urban","forest"),times=c(8,8,8))
IL6data<-data.frame(IL6,IL6groups)
s9<- ggplot(IL6data, aes(factor(IL6groups),IL6))
s9<-s9+geom_boxplot(aes(fill=factor(IL6groups)),show.legend = F)+ylab("IL-6 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s9
kruskal.test(IL6~IL6groups,IL6data)

IL66<-c(2.23,
        1.44,
        3.10,
        1.99,
        2.15,
        2.15,
        22.37,
        21.12,
        15.05,
        2.70,
        9.33,
        13.41,
        1.99,
        1.36,
        2.15,
        0.09,
        0.73,
        1.83,
        1.04)
IL66groups=rep(c("urban","forest","no soil"),times=c(8,8,3))
IL66data<-data.frame(IL66,IL66groups)
s69<- ggplot(IL66data, aes(factor(IL66groups),IL66))
s69<-s69+geom_boxplot(aes(fill=factor(IL66groups)),show.legend = F)+ylab("IL-6 (pg/mL)")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
s69
kruskal.test(IL66~IL66groups,IL66data)

##wrangle
colnames(IL6data)<-c("cytokine", "treatment")
colnames(IL66data)<-c("cytokine", "treatment")
IL6data$age<-3
IL66data$age<-6
IL6Total <- rbind(IL6data,IL66data)

##plot 
IL6Totalplot<- ggplot(IL6Total, aes(factor(treatment),cytokine)) + 
  geom_boxplot(aes(fill=factor(age)), show.legend=T)+ylab("serum cytokine concentration (pg/mL)")+
  xlab("")+
  guides(fill=guide_legend(title="Age(weeks)"))+
  geom_jitter()+
  ggtitle("IL-6")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
IL6Totalplot

##subset
IL6F<-subset(IL6Total, treatment=="forest")
IL6U<-subset(IL6Total, treatment=="urban")
IL6N<-subset(IL6Total, treatment=="no soil")

##wilcoxon, if p-value is less than .05 significance level, reject null hypothesis 
wilcox.test(cytokine ~ age, data=IL6F) ## 0.43
wilcox.test(cytokine ~ age, data=IL6U) ## 0.16
wilcox.test(cytokine ~ age, data=IL6N) ## 0.01 reject null!

##multiplot function
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
##make multiplot
multiplot(EotaxinTotalplot,GCSFTotalplot,IL1aTotalplot,IL1bTotalplot,IL10Totalplot,MCP1Totalplot,TnfaTotalplot,IL13Totalplot,IL6Totalplot,cols = 3)
##writeout
pdffile <- "EVE_Dec10.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
EotaxinTotalplot +theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
GCSFTotalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
IL1aTotalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
IL1bTotalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
IL10Totalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
MCP1Totalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
TnfaTotalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
IL13Totalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
IL6Totalplot+ theme(text=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y=element_text(size=14))
dev.off()

MPO<-c(6,2,1,1,0,0,0,0,0,0,4,4,3,1,0,0,0,0,0,0,0,1,0)
MPOgroups=rep(c("forest","urban","no soil"),times=c(10,8,5))
MPOdata<-data.frame(MPO,MPOgroups)
b1<- ggplot(MPOdata, aes(factor(MPOgroups),MPO))
b1<-b1+geom_boxplot(aes(fill=factor(MPOgroups)),show.legend = F)+ylab("MPO+ cell counts/colon section")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
b1
kruskal.test(MPO~MPOgroups,MPOdata)

MPO6<-c(6,0,0,0,0,0,5,2,0,0,0,1,0,0)
MPO6groups=rep(c("forest","urban","no soil"),times=c(6,5,3))
MPO6data<-data.frame(MPO6,MPO6groups)
b2<- ggplot(MPO6data, aes(factor(MPO6groups),MPO6))
b2<-b2+geom_boxplot(aes(fill=factor(MPO6groups)),show.legend = F)+ylab("MPO+ cell counts/colon section")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
b2
kruskal.test(MPO6~MPO6groups,MPO6data)

CD8<-c(5,0,0,2,1,3,0,0,10,5,3,1,2,0,1,0,2,3,3,4,2,0,1,0)
CD8groups=rep(c("forest","urban","no soil"),times=c(8,8,8))
CD8data<-data.frame(CD8,CD8groups)
b3<- ggplot(CD8data, aes(factor(CD8groups),CD8))
b3<-b3+geom_boxplot(aes(fill=factor(CD8groups)),show.legend = F)+ylab("CD8+ T cell counts/colon section")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
b3
kruskal.test(CD8~CD8groups,CD8data)

CD86<-c(2,0,1,0,3,0,5,2,1,0,3,1,0,0)
CD86groups=rep(c("forest","urban","no soil"),times=c(6,5,3))
CD86data<-data.frame(CD86,CD86groups)
b4<- ggplot(CD86data, aes(factor(CD86groups),CD86))
b4<-b4+geom_boxplot(aes(fill=factor(CD86groups)),show.legend = F)+ylab("CD8+ T cell counts/colon section")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
b4
kruskal.test(CD86~CD86groups,CD86data)

F480<-c(1,5,1,1,1,5,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
F480groups=rep(c("forest","urban","no soil"),times=c(11,10,5))
F480data<-data.frame(F480,F480groups)
b5<- ggplot(F480data, aes(factor(F480groups),F480))
b5<-b5+geom_boxplot(aes(fill=factor(F480groups)),show.legend = F)+ylab("F4/80+ cell counts/colon section")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
b5
kruskal.test(F480~F480groups,F480data)

F4806<-c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0)
F4806groups=rep(c("urban","forest","no soil"),times=c(7,12,3))
F4806data<-data.frame(F4806,F4806groups)
b6<- ggplot(F4806data, aes(factor(F4806groups),F4806))
b6<-b6+geom_boxplot(aes(fill=factor(F4806groups)),show.legend = F)+ylab("F4/80+ cell counts/colon section")+xlab("")+ geom_jitter()+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
b6
kruskal.test(F4806~F4806groups,F4806data)


##output file
pdffile <- "qPCR.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
multiplot(r,r2,r3,r4,r5,r6,cols = 3)
dev.off()

pdffile <- "qPCR_EVE.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
multiplot(s,s2,s3,s4,s5,s6,s7,s8,s9,cols = 3)
dev.off()

pdffile <- "qPCR_EVE_6weeks.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
multiplot(s61,s62,s63,s64,s65,s66,s67,s68,s69,cols = 3)
dev.off()

pdffile <- "MPO_counts.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=5, height=3)
multiplot(b1,b2,cols = 2)
dev.off()

pdffile <- "CD8_counts.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=5, height=3)
multiplot(b3,b4,cols = 2)
dev.off()

pdffile <- "F480_counts.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=5, height=3)
multiplot(b5,b6,cols = 2)
dev.off()

