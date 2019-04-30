##GC data total visualized with grouped barplot
##package ggplot2, PMCMR

##here I'm making my dataframe from scratch, but you can upload a csv and format it alternatively

##3 weeks
SCFAtotal3<-c(0.09914378,
              0.14973298,
              0.07354353,
              0.16855557,
              0.08054611,
              0.14443624,
              0.05936801,
              0.23114945,
              0.14175629,
              0.09077993,
              0.02384603,
              0.36585470,
              0.11754847,
              0.09826760,
              0.07864670,
              0.05936801,
              0.15324654,
              0.33452077,
              0.37983918,
              0.16460105,
              0.109044036,
              0.129223997,
              0.15324654,
              0.33452077,
              0.37983918,
              0.09077993,
              0.10301167,
              0.16160004,
              0.04145239,
              0.50567814,
              0.13408589,
              0.10301167,
              0.161600039,
              0.167613656)
groups_total3=rep(c("no soil","urban","forest"),times=c(10,12,12))
total3data<-data.frame(SCFAtotal3,groups_total3)

total3plot<- ggplot(total3data, aes(factor(groups_total3),SCFAtotal3))
total3plot<-total3plot+geom_boxplot(aes(fill=factor(groups_total3)),show.legend = F)+ylab("Total cecal SCFA (% weight)")+xlab("")+ggtitle("3 weeks")+ geom_jitter()
total3plot<-total3plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
total3plot

kruskal.test(SCFAtotal3~groups_total3,total3data)
attach(total3data)
dunnA<-posthoc.kruskal.dunn.test(x=total3data$SCFAtotal3, g=total3data$groups_total3, p.adjust.methods="bonferroni")
summary(dunnA)

#6 weeks
SCFAtotal6<-c(0.04145239,
              0.50567814,
              0.23114945,
              0.16761366,
              0.13408589,
              0.14175629,
              0.16855557,
              0.08054611,
              0.14443624,
              0.02384603,
              0.36585470,
              0.07864670,
              0.05680572,
              0.10290811,
              0.16460105,
              0.08497672,
              0.10904404,
              0.01338279,
              0.09914378,
              0.14973298,
              0.07354353,
              0.11754847,
              0.09826760,
              0.43778832,
              0.12922400,
              0.43778832,
              0.05680572,
              0.10290811,
              0.08497672,
              0.01338279)
groups_total6=rep(c("urban","forest","no soil"),times=c(12,12,6))
total6data<-data.frame(SCFAtotal6,groups_total6)

total6plot<- ggplot(total6data, aes(factor(groups_total6),SCFAtotal6))
total6plot<-total6plot+geom_boxplot(aes(fill=factor(groups_total6)),show.legend = F)+ylab("Total cecal SCFA (% weight)")+xlab("")+ggtitle("6 weeks")+ geom_jitter()
total6plot<-total6plot+theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
total6plot

kruskal.test(SCFAtotal6~groups_total6,total6data)
attach(total6data)
dunnA<-posthoc.kruskal.dunn.test(x=total6data$SCFAtotal6, g=total6data$groups_total6, p.adjust.methods="bonferroni")
summary(dunnA)

pdffile <- "GC_total.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=6)
multiplot(total3plot,total6plot,cols=2)
dev.off()

multiplot(total3plot,total6plot,cols=2)
