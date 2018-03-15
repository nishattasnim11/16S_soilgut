##GC data

##3 weeks
acetic3<-c(0.062247863,
           0.091791296,
           0.040214006,
           0.095371687,
           0.048701625,
           0.088554799,
           0.058522279,
           0.184519795,
           0.029909198,
           0.046000445,
           0.01285672,
           0.205857592,
           0.075794087,
           0.051785002,
           0.038839371,
           0.037436093,
           0.14246034,
           0.03661173,
           0.194558298,
           0.013015629,
           0.090674193,
           0.079445456,
           0.104356367,
           0.195424148,
           0.220211684,
           0.052914822,
           0.053935689,
           0.089291143,
           0.09898163,
           0.06903061,
           0.063006858,
           0.1432067,
           0.121507794,
           0.102563537)
acetic3groups=rep(c("no soil","urban","forest"),times=c(10,12,12))
acetic3data<-data.frame(acetic3,acetic3groups)
acetic3plot<- ggplot(acetic3data, aes(factor(acetic3groups),acetic3))
acetic3plot<-acetic3plot+geom_boxplot(aes(fill=factor(acetic3groups)),show.legend = F)+ylab("% weight")+xlab("")+ggtitle("Acetic Acid")+ geom_jitter()
acetic3plot<-acetic3plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))
kruskal.test(acetic3~acetic3groups,acetic3data)
attach(acetic3data)
dunnA<-posthoc.kruskal.dunn.test(x=acetic3data$acetic3, g=acetic3data$acetic3groups, p.adjust.methods="bonferroni")
summary(dunnA)

propionic3<-c(0.0099414,
              0.010891286,
              0.008142598,
              0.013050903,
              0.012429083,
              0.010753509,
              0.004378746,
              0.024510193,
              0.003389766,
              0.003656354,
              0.004757044,
              0.030564658,
              0.010312456,
              0.009622902,
              0.010418916,
              0.005284216,
              0.014903612,
              0.003741835,
              0.01771716,
              0.002519542,
              0.009493854,
              0.008884289,
              0.016038638,
              0.02319639,
              0.024986103,
              0.009596128,
              0.011921697,
              0.011195844,
              0.014629573,
              0.00521355,
              0.006839255,
              0.020657433,
              0.016207873,
              0.014522346)
propionic3groups=rep(c("no soil","urban","forest"),times=c(10,12,12))
propionic3data<-data.frame(propionic3,propionic3groups)
propionic3plot<- ggplot(propionic3data, aes(factor(propionic3groups),propionic3))
propionic3plot<-propionic3plot+geom_boxplot(aes(fill=factor(propionic3groups)),show.legend = F)+ylab("% weight")+xlab("")+ggtitle("Propionic Acid")+ geom_jitter()
propionic3plot<-propionic3plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))

kruskal.test(propionic3~propionic3groups,propionic3data)
##need to attach for downstream steps
attach(propionic3data)
##dunns pairwise tests for Observed (note: best to use p.adjust.method= “bonferroni” or sth)
dunnP<-posthoc.kruskal.dunn.test(x=propionic3data$propionic3, g=propionic3data$propionic3groups, p.adjust.methods="bonferroni")
#check results
summary(dunnP)

butyric3<-c(0.026954512,
            0.047050395,
            0.025186929,
            0.060132977,
            0.019415403,
            0.045127936,
            0.035064101,
            0.04649454,
            0.009751397,
            0.02050265,
            0.006232264,
            0.129432447,
            0.031441927,
            0.036859696,
            0.029388408,
            0.016647698,
            0.072828215,
            0.015597832,
            0.1128756,
            0.004633845,
            0.028142466,
            0.038806492,
            0.032851537,
            0.11590023,
            0.134641393,
            0.028268984,
            0.037154284,
            0.061113052,
            0.045711007,
            0.028202721,
            0.035673404,
            0.08048414,
            0.067115396,
            0.05899724)
butyric3groups=rep(c("no soil","urban","forest"),times=c(10,12,12))
butyric3data<-data.frame(butyric3,butyric3groups)
butyric3plot<- ggplot(butyric3data, aes(factor(butyric3groups),butyric3))
butyric3plot<-butyric3plot+geom_boxplot(aes(fill=factor(butyric3groups)),show.legend = F)+ylab("% weight")+xlab("")+ggtitle("Butyric Acid")+ geom_jitter()
butyric3plot<-butyric3plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))

kruskal.test(butyric3~butyric3groups,butyric3data)
##need to attach for downstream steps
attach(butyric3data)
##dunns pairwise tests for Observed (note: best to use p.adjust.method= “bonferroni” or sth)
dunnB<-posthoc.kruskal.dunn.test(x=butyric3data$butyric3, g=butyric3data$butyric3groups, p.adjust.methods="bonferroni")
#check results
summary(dunnB)

##6 weeks
acetic6<-c(0.02342292,
           0.268520002,
           0.119314055,
           0.090831149,
           0.082491658,
           0.08489663,
           0.12836715,
           0.121077521,
           0.077325433,
           0.194225493,
           0.210507242,
           0.201659386,
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
           0.13160898,
           0.071275487,
           0.274593971,
           0.096723106,
           0.035117897,
           0.052347315,
           0.052394675)
acetic6groups=rep(c("urban","forest","no soil"),times=c(12,12,6))
acetic6data<-data.frame(acetic6,acetic6groups)
acetic6plot<- ggplot(acetic6data, aes(factor(acetic6groups),acetic6))
acetic6plot<-acetic6plot+geom_boxplot(aes(fill=factor(acetic6groups)),show.legend = F)+ylab("% weight")+xlab("")+ggtitle("Acetic Acid")+geom_jitter()
acetic6plot<-acetic6plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))

kruskal.test(acetic6~acetic6groups,acetic6data)
attach(acetic6data)
dunnA6<-posthoc.kruskal.dunn.test(x=acetic6data$acetic6, g=acetic6data$acetic6groups, p.adjust.methods="bonferroni")
summary(dunnA6)

propionic6<-c(0.006411566,
              0.040816552,
              0.023102838,
              0.020636736,
              0.014653116,
              0.015710545,
              0.024772991,
              0.013458956,
              0.017301995,
              0.035195763,
              0.028530234,
              0.031307162,
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
              0.014843349,
              0.015400307,
              0.038285025,
              0.009813442,
              0.003596463,
              0.005750894,
              0.004620294)
propionic6groups=rep(c("urban","forest","no soil"),times=c(12,12,6))
propionic6data<-data.frame(propionic6,propionic6groups)
propionic6plot<- ggplot(propionic6data, aes(factor(propionic6groups),propionic6))
propionic6plot<-propionic6plot+geom_boxplot(aes(fill=factor(propionic6groups)),show.legend = F)+ylab("% weight")+xlab("")+ggtitle("Propionic Acid")+ geom_jitter()
propionic6plot<-propionic6plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))

kruskal.test(propionic6~propionic6groups,propionic6data)
attach(propionic6data)
dunnP6<-posthoc.kruskal.dunn.test(x=propionic6data$propionic6, g=propionic6data$propionic6groups, p.adjust.methods="bonferroni")
summary(dunnP6)

butyric6<-c(0.011617904,
            0.196341584,
            0.08873256,
            0.056145771,
            0.036941114,
            0.041149117,
            0.080650054,
            0.057707458,
            0.032319282,
            0.0790243,
            0.101457921,
            0.10515647,
            0.018277685,
            0.02931816,
            0.05401128,
            0.021838374,
            0.029633223,
            0.003357772,
            0.064938254,
            0.121031426,
            0.00200193,
            0.10639277,
            0.007736106,
            0.041238869,
            0.042548204,
            0.124909324,
            0.051441451,
            0.018649916,
            0.026222813,
            0.024215657)
butyric6groups=rep(c("urban","forest","no soil"),times=c(12,12,6))
butyric6data<-data.frame(butyric6,butyric6groups)
butyric6plot<- ggplot(butyric6data, aes(factor(butyric6groups),butyric6))
butyric6plot<-butyric6plot+geom_boxplot(aes(fill=factor(butyric6groups)),show.legend = F)+ylab("% weight")+xlab("")+ggtitle("Butyric Acid")+ geom_jitter()
butyric6plot<-butyric6plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=16))

kruskal.test(butyric6~butyric6groups,butyric6data)
attach(butyric6data)
dunnB6<-posthoc.kruskal.dunn.test(x=butyric6data$butyric6, g=butyric6data$butyric6groups, p.adjust.methods="bonferroni")
summary(dunnB6)

pdffile <- "GC_Nov19.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
multiplot(acetic3plot,propionic3plot,butyric3plot,acetic6plot,propionic6plot,butyric6plot,cols = 2)
dev.off()



