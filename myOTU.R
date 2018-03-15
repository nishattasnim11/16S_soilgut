##Here I will output a csv file of my OTU table with taxonomy and metadata fields

#let's load the packages first 
.load_packages <- c("phyloseq", "ggplot2", "plyr", "vegan","DESeq2") 
sapply(c(.load_packages), require, character.only = TRUE)

#I'm setting my working directory to Desktop setwd("~/Desktop")
myOTU<-import_biom(BIOMfilename = "/Users/nishattasnim/Desktop/Dec 2017/May 30th 2017/otu_943cutoff_with_sample_metadata_json.biom", treefilename = "/Users/nishattasnim/Desktop/May 30th 2017/rep_set.tre",parseFunction = parse_taxonomy_greengenes)

#check the phyloseq object myOTU, check: https://joey711.github.io/phyloseq/preprocess.html
#rename the columns of the tax table
colnames(tax_table(myOTU)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
##turn group variable from sample data into factor 
factor(x = (sample_data(myOTU)$group), levels = c("forest", "urban", "no soil"),exclude = "na")
#can check rank names with rank_names(myOTU)

wh0= genefilter_sample(myOTU, filterfun_sample(function(x) x > 5), A=0.5*nsamples(myOTU))
myOTU1= prune_taxa(wh0, myOTU)
myOTU1= transform_sample_counts(myOTU1, function(x) 1E6 * x/sum(x))
phylum.sum = tapply(taxa_sums(myOTU1), tax_table(myOTU1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
myOTU1= prune_taxa((tax_table(myOTU1)[, "Phylum"] %in% top5phyla), myOTU1)

myOTU1_gut=subset_samples(myOTU1, env_matter=="ENVO:intestine environment")
myOTU1_gut_3weeks=subset_samples(myOTU1_gut, age_weeks=="3")
myOTU1_gut_6weeks=subset_samples(myOTU1_gut, age_weeks=="6")

p = plot_bar(gutOTU, "Phylum", fill="Phylum", facet_grid=group~age_weeks)
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

##Beta plots
myOTU1_gut_3weeks.ord <- ordinate(myOTU1_gut_3weeks, "NMDS", "bray")
ord3weeeks=plot_ordination(myOTU1_gut_3weeks,myOTU1_gut_3weeks.ord,color="group")+stat_ellipse()
ord3weeeks=ord3weeeks+geom_point(aes(color = group), alpha = 0.7, size = 4)
ord3weeeks=ord3weeeks+theme(plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=10) , legend.text=element_text(size=9))

print(ord3weeeks)

stressplot(myOTU1_gut_3weeks.ord)

myOTU1_gut_6weeks.ord <- ordinate(myOTU1_gut_6weeks, "NMDS", "bray")
ord6weeeks=plot_ordination(myOTU1_gut_6weeks, myOTU1_gut_6weeks.ord, color="group")+stat_ellipse()
ord6weeeks=ord6weeeks+geom_point(aes(color = group), alpha = 0.7, size = 5)
ord6weeeks=ord6weeeks+theme(plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=10) , legend.text=element_text(size=9))

print(ord6weeeks)
stressplot(myOTU1_gut_6weeks.ord)

#Next, perform taxa-wise filtering: remove taxa not seen more than 3 times in at least 20% of the samples 
myOTU_taxFiltered = filter_taxa(myOTU, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
myOTU_taxFiltered

#check filtered table myOTU_taxaFiltered (now has 2457 taxa of the original ~14092 taxa)
# I plotted this heat map for fun: plot_heatmap(myOTU_taxFiltered)

#finally, let's write out the file so we can send it to Bod
dat <- psmelt(myOTU_taxFiltered)
write.csv(dat, file='otus-with-sample-data.csv')

###making graphs for thesis document from this session

soilOTU = subset_samples(myOTU, env_matter=="ENVO:soil")
gutOTU = subset_samples(myOTU, env_matter=="ENVO:intestine environment")

gutOTU=prune_species(speciesSums(gutOTU) > 0, gutOTU)
soilOTU=prune_species(speciesSums(soilOTU) > 0, soilOTU)

gutOTU_3weeks=subset_samples(gutOTU, age_weeks=="3")
(sample_data(gutOTU_3weeks)$group)=ordered(x = (sample_data(gutOTU_3weeks)$group), levels = c("forest", "urban", "no soil"))

gutOTU_6weeks=subset_samples(gutOTU, age_weeks=="6")

gutOTU_3weeks_forest=subset_samples(gutOTU_3weeks, group="forest")
gutOTU_3weeks_urban=subset_samples(gutOTU_3weeks, group="urban")
gutOTU_3weeks_nosoil=subset_samples(gutOTU_3weeks, group="no soil")

gutOTU_6weeks_forest=subset_samples(gutOTU_6weeks, group="forest")
gutOTU_6weeks_urban=subset_samples(gutOTU_6weeks, group="urban")
gutOTU_6weeks_nosoil=subset_samples(gutOTU_6weeks, group="no soil")

myOTU_taxFiltered_soil=subset_samples(myOTU_taxFiltered, env_matter=="ENVO:soil")
myOTU_taxFiltered_gut=subset_samples(myOTU_taxFiltered, env_matter=="ENVO:intestine environment")
myOTU_taxFiltered_gut_3weeks=subset_samples(myOTU_taxFiltered_gut, age_weeks=="3")
myOTU_taxFiltered_gut_6weeks=subset_samples(myOTU_taxFiltered_gut, age_weeks=="6")

taxasoilplot2=plot_taxa_summary(myOTU_taxFiltered_soil, "Phylum")
taxasoilplot2=taxasoilplot2+theme(plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=10) , legend.text=element_text(size=9))
taxasoilplot2

taxagutplot2_3weeks=plot_taxa_summary(myOTU_taxFiltered_gut_3weeks, "Phylum")
taxagutplot2_3weeks=taxagutplot2_3weeks+theme(plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=10) , legend.text=element_text(size=9))
taxagutplot2_3weeks


taxagutplot2_6weeks=plot_taxa_summary(myOTU_taxFiltered_gut_6weeks, "Phylum")
taxagutplot2_6weeks=taxagutplot2_6weeks+theme(plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=10) , legend.text=element_text(size=9))
taxagutplot2_6weeks


grid.arrange(taxagutplot2_3weeks, taxagutplot2_6weeks, ncol = 2)

library(grid)
library(gridExtra)
pdf("taxa_summary_Phylum2_gut_age.pdf",width = 7,height = 5)
grid.arrange(taxagutplot2_3weeks, taxagutplot2_6weeks, ncol = 2)
dev.off()

##richness plot 
plot_bar(myOTU_taxFiltered_soil, fill="Phylum")+ theme(plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=10) , legend.text=element_text(size=9))

plot_heatmap(myOTU_taxFiltered_gut_6weeks)

###Let's try some differential abundance with DESeq2
diagddsGut = phyloseq_to_deseq2(gutOTU, ~ group)
diagddsGut = DESeq(diagddsGut, test="Wald", fitType="parametric")

pdffile <- "abundance_3&6weeks.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
par(mfrow=c(2,1))
plot_bar(myOTU_taxFiltered_gut_3weeks, "group", fill="Genus", facet_grid=~Phylum)
plot_bar(myOTU_taxFiltered_gut_6weeks, "group", fill="Genus", facet_grid=~Phylum)
dev.off()

plot_bar(myOTU_taxFiltered_gut_3weeks, "", fill="Species", facet_grid=~Genus)

GutFirm3=subset_taxa(myOTU_taxFiltered_gut_3weeks,Phylum=="Firmicutes")
plot_bar(GutFirm3, "group", fill="Family", facet_grid=~Genus)
GutFirm3_Lacto=subset_taxa(myOTU_taxFiltered_gut_3weeks,Genus=="Lactobacillus")
plot_richness(GutFirm3, "group") + geom_boxplot()
plot_bar(GutFirm3, "group")

g1=plot_richness(gutOTU_3weeks,x = "group",color = "group",measures = c("Observed","Chao1","Shannon","Simpson","InvSimpson"))
g1 + geom_point(size=5, alpha=0.7)
g1$layers
g1$layers <- g1$layers[-1]
g1 + geom_point(size=5, alpha=0.7)

##Richness plots 
pdffile <- "alpha_3&6weeks_observed.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
par(mfrow=c(2,1))
plot_richness(gutOTU_3weeks,x = "group",color ="group",measures = c("Observed"))+geom_boxplot(alpha=0.7)+theme(text=element_text(size=16), axis.text.x = element_text(angle = 360, hjust=0.5))
plot_richness(gutOTU_6weeks,x = "group",color ="group",measures = c("Observed"))+geom_boxplot(alpha=0.7)+theme(text=element_text(size=16), axis.text.x = element_text(angle = 360, hjust=0.5))
dev.off()

##prune out OTUs that occur 0 times 
#GPalphaOTU <- prune_taxa(taxa_sums(alphaOTU) > 0, alphaOTU)
##obtain estimates
alphaDiversity3<-estimate_richness(gutOTU_3weeks,measures = c("Observed","Chao1","Shannon","Simpson","InvSimpson"))
alphaDiversity6<-estimate_richness(gutOTU_6weeks,measures = c("Observed","Chao1","Shannon","Simpson","InvSimpson"))

##add estimates to sample data of phyloseq object
data3 <- cbind(sample_data(gutOTU_3weeks), alphaDiversity3)
data6 <- cbind(sample_data(gutOTU_6weeks), alphaDiversity6)

##make category of analysis a factor
data3$group <- as.factor(data3$group)
data6$group <- as.factor(data6$group)

##do kruskal test on group category (use whichever metric you like, here I used Observed richness)
kruskal.test(Observed ~ group, data=data3)
kruskal.test(Chao1 ~ group, data=data3)
kruskal.test(Shannon ~ group, data=data3)
kruskal.test(Simpson ~ group, data=data3)
kruskal.test(InvSimpson ~ group, data=data3)

kruskal.test(Observed ~ group, data=data6)
kruskal.test(Chao1 ~ group, data=data6)
kruskal.test(Shannon ~ group, data=data6)
kruskal.test(Simpson ~ group, data=data6)
kruskal.test(InvSimpson ~ group, data=data6)

##need to attach for downstream steps
attach(data3)
##dunns pairwise tests for Observed (note: best to use p.adjust.method= “bonferroni” or sth)
dunn1<-posthoc.kruskal.dunn.test(x=Observed, g=group, p.adjust.methods="bonferroni")
#check results
summary(dunn1)

##beta diversity plots
pdffile <- "beta_3&6weeks.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
par(mfrow=c(2,1))
print(ord3weeeks)+theme(text=element_text(size=16))
print(ord6weeeks)+theme(text=element_text(size=16))
dev.off()

library(magrittr)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
###link: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
###Stacked barplot
gutOTU_StackedBarplot= merge_phyloseq(gutOTU_3weeks,gutOTU_6weeks)
gutOTU_phylum <- gutOTU_StackedBarplot %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


gutOTU_phylumRel <- gutOTU_StackedBarplot %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)/nsamples(gutOTU_StackedBarplot)} )  %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

sample_data(gutOTU_phylum)$group <- as.factor(x = (sample_data(gutOTU_phylum)$group))

Barplot_gut<-ggplot(gutOTU_phylum, aes(x = group, y = Abundance, fill = Phylum)) + 
  facet_grid(age_weeks~.) +
  geom_bar(stat = "identity",position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = c("forest", "urban", "no soil"),
    labels = c("Forest", "Urban", "No soil"), 
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
print(Barplot_gut)


levels(sample_data(gutOTU_phylum)$group)
##set levels to re-do grouping of barplot

sample_data(mydata)$BMIClass = factor(sample_data(mydata)$BMIClass, levels = c("SeverelyUnderweight","Underweight","Normal","Overweight", "Obese"))
##permanova
set.seed(1)
gut3_bray <- phyloseq::distance(myOTU1_gut_3weeks, method = "bray")
gut3sampledf<-data.frame(sample_data(myOTU1_gut_3weeks))
adonis(gut3_bray ~ group, data = gut3sampledf)
anosim(gut3_bray ~ group, data = gut3sampledf)
anosim(phyloseq::distance(myOTU1_gut_3weeks, method = "bray"),sample_data(myOTU1_gut_3weeks)$group)

gut6_bray <- phyloseq::distance(myOTU1_gut_6weeks, method = "bray")
gut6sampledf<-data.frame(sample_data(myOTU1_gut_6weeks))
adonis(gut6_bray ~ group, data = gut6sampledf)
anosim(phyloseq::distance(myOTU1_gut_6weeks, method = "bray"),sample_data(myOTU1_gut_6weeks)$group)


# Homogeneity of dispersion test
beta3 <- betadisper(gut3_bray, gut3sampledf$group)
permutest(beta3)
boxplot(beta3)
TukeyHSD(beta3)

beta6 <- betadisper(gut6_bray, gut6sampledf$group)
permutest(beta6)
boxplot(beta6)
TukeyHSD(beta6)

##Now for soil composition
soilOTU_phylum <- soilOTU %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)  

Barplot_soil<-ggplot(soilOTU_phylum, aes(x = ph, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity",position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = c("7.35", "5.27"),
    labels = c("A&W","Chute Lake"), 
    drop = FALSE
  )+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n")+
  theme(text=element_text(size=12), axis.text.x= element_text(size=12), axis.text.y=element_text(size=12))
print(Barplot_soil)

factor(x = (sample_data(gutOTU_3weeks)$group), levels = c("forest", "urban", "no soil"))


pdffile <- "phylumcomposition_gut.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
print(Barplot_gut)
dev.off()

pdffile <- "phylumcomposition_soil.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
print(Barplot_soil)
dev.off()

pdffile <- "beta_3&6weeks_ellipse2.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
par(mfrow=c(2,1))
print(ord3weeeks)
print(ord6weeeks)
dev.off()

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(myOTU))


# Histogram of sample read counts
pdffile <- "SampleSequencingDepth.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
dev.off()  

smin <- min(sample_sums(myOTU))
smean <- mean(sample_sums(myOTU))
smax <- max(sample_sums(myOTU))

###multiplot function
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

##more relative abundance plots
library("phyloseq")
library("data.table")
library("ggplot2")

fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

plot_taxa_summary = function(physeq, Rank, GroupBy = NULL){
  # Get taxa summary table 
  dt1 = summarize_taxa(physeq, Rank = Rank, GroupBy = GroupBy)
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]], 
                          levels = rev(dt1[[Rank]]))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]
  
  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5)
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin))
  }
  return(pRank)
}

# Test
data("gutOTU")
plot_taxa_summary(gutOTU_3weeks_forest, "Phylum")
gutOTU_3weeks_forest=subset_samples(gutOTU_3weeks, group="forest")
gutOTU_3weeks_urban=subset_samples(gutOTU_3weeks, group="urban")
gutOTU_3weeks_nosoil=subset_samples(gutOTU_3weeks, group="no soil")

plot_taxa_summary(gutOTU, "Phylum", GroupBy = "group")
summarize_taxa(gutOTU, "Phylum")
summarize_taxa(gutOTU, "Phylum", "age_weeks")

pdffile <- "taxa_summaries_3&6weeks.pdf"
pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
par(mfrow=c(3,2))
plot_taxa_summary(gutOTU_3weeks_forest, "Phylum")
plot_taxa_summary(gutOTU_3weeks_urban, "Phylum")
plot_taxa_summary(gutOTU_3weeks_nosoil, "Phylum")
plot_taxa_summary(gutOTU_6weeks_forest, "Phylum")
plot_taxa_summary(gutOTU_6weeks_urban, "Phylum")
plot_taxa_summary(gutOTU_6weeks_nosoil, "Phylum")
dev.off()

pdf(paste(pdffile,sep=""), paper="letter", width=8, height=8)
plot_richness(soilOTU)
dev.off()
