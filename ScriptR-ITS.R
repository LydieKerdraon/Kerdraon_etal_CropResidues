#this script can be use on ITS data 
#Data used for these analyses are given in the R environment Data_ITS.RData

####################Environnement
# Haplo : Haplotypes list with 298 samples. (149 biological samples* 2 technical replicas)
# newname : file used to rename sample names and delete haplotypes not present in the two technical replicas ("Uncertain presences")
# sample.dataITS : description of biological samples
# tax_tableITS: taxonomic assignement of haplotypes performed on Unite 7.1 database. In this script, we include the taxonomy table made with : taxaITS <- assignTaxonomy(haplo.tableITSf, "path_to_database", multithread=TRUE)
# Temp1 & Temp2 : haplotypes presence tables. Files used for decomposition of beta-diversity with betapart Package

######## Packages required + custom functions
Packages <- c("phyloseq", "data.table", "ggplot2", "plyr","dplyr","reshape2","grid",
              "gridExtra","scales","dplyr", "ggpubr","vegan","multcompView","rcompanion","betapart")

lapply(Packages, library, character.only = TRUE)


####### "Uncertain presence" removing

Haplo<-Haplo[order(match(rownames(newname), Haplo[1])), ]
Haplo[1]<-newname[1]

Haplo[Haplo==0]<-NA
NewHaplo.Table<-aggregate(. ~ X,Haplo, sum, na.action = na.pass)
NewHaplo.Table[is.na(NewHaplo.Table)] <- 0

NewHaplo.Table.with.rownames <- data.frame(NewHaplo.Table[,-1], row.names=NewHaplo.Table[,1])
NewHaplo.Table2<-as.matrix(NewHaplo.Table.with.rownames,row.names=TRUE, header=T)

haplo.tableITS<-t(NewHaplo.Table2)
haplo.tableITSf = haplo.tableITS[ rowSums(haplo.tableITS)!=0, ] 


##############
######## Phyloseq object preparation [phyITS]
tax.tableITS <- replace(taxaITS, is.na(taxaITS), "Unclassified")
sample.dataITS$Month <- factor(sample.dataITS$Month, levels = c("Jul. ", "Oct. ", "Dec. ", "Feb. ","May"))

phyITS <- phyloseq(otu_table(haplo.tableITSf, taxa_are_rows = TRUE),
                 tax_table(tax.tableITS),
                 sample_data(sample.dataITS))

phyITS

#Removing haplo=0
condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phyITS, condition, 1)
phyITS<-prune_taxa(taxaToKeep, phyITS)
phyITS

#Removing unclassified haplotypes / Anthophyta, Arthropoda, Cercozoa,unclassified_Plantae
table(tax_table(phyITS)[, "Phylum"])
phyITS <- subset_taxa(phyITS, !Phylum %in% c("", "p__Anthophyta", "p__Arthropoda", "p__Cercozoa", "p__unclassified_Plantae", "Unclassified"))
phyITS

######## Reads number and normalization
TotalReads = sample_sums(phyITS)

#sequencing depth
SeqDepth = data.table(as(sample_data(phyITS), "data.frame"),
                  TotalReads = sample_sums(phyITS), keep.rownames = TRUE)
pSeqDepth = ggplot(SeqDepth, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
print(pSeqDepth)

#Normalization
phyITSnorm <- transform_sample_counts(phyITS, function(x) round(x/sum(x) *1000000 ))
#write.csv(cbind((otu_table(phyITSnorm)),tax_table(phyITSnorm)), "ITS-data.csv")


######## Alpha-diversity analyses : Shannon index
cols <- c("Oilseed Rape" = "#ea0b16", "Wheat (monoculture)" = "#29ce37", "Wheat (Rotation)" = "#007bcf")
measures=c("Shannon")
	#various possible indexes : "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"

RichPlotITS<-plot_richness(phyITSnorm, x="Month", measures = measures)+ 
	stat_compare_means(label.y=c(0.2),label.x.npc = 0.2, method="kruskal.test")+
	geom_boxplot(outlier.shape =NA, aes(fill=Plant.Succession), lwd=0.15, width = 0.3) +#boxplot
	theme_gray(base_size = 10) + ylab(paste(measures,"index (Fungal diversity)")) + xlab("")

RichPlotITS+facet_grid(Growing.Season~Plant.Succession)


##----Wilcoxon rank sum test on Shannon index -- Example for 1st growing season
richnessITS<-estimate_richness(phyITSnorm,measures = measures)

test<-sample.dataITS[match(rownames(richnessITS), rownames(sample.dataITS)), ]
compareRichITS<-cbind(test,richnessITS)


Year1<-subset(compareRichITS, Growing.Season=="2015-2016")

OilseedITS1 <- subset(Year1, Plant.Succession=="Oilseed Rape") 
WheatMITS1 <- subset(Year1, Plant.Succession=="Wheat (monoculture)")
WheatRITS1 <- subset(Year1, Plant.Succession=="Wheat (Rotation)")

#example : wilcoxon test between plant succession
compare_means(data=Year1,Shannon~Plant.Succession,method = "wilcox.test")

#example : pairwise wilcoxon for wheat monoculture, 1st growing season
WheatM1.ITS.pairwise<-pairwise.wilcox.test(WheatMITS1$Shannon,WheatMITS1$Month, p.adjust.method = "none")
WheatM1.ITS.pvalue = fullPTable(WheatM1.ITS.pairwise$p.value)
multcompLetters(WheatM1.ITS.pvalue) #require : rcompanion package

######## Beta-diversity analyses : Multidimensional scaling on Bray-curtis index
OTU.ordITS<-ordinate(phyITSnorm,"MDS","bray")
MdsITS=plot_ordination(phyITSnorm,OTU.ordITS,type="sample",color="Plant.Succession",shape="Growing.Season")+
	scale_shape_manual(values = c(23,16,17)) + geom_point(size=3)+
	facet_grid(.~Month)+ 
	scale_color_manual(values=cols)+
	theme(panel.spacing = unit(1, "lines"))+
	theme_bw()
MdsITS
##----Permanova test on Bray-curtis index -- Example for "Plant" factor on overall dataset
adonis(distance(phyITSnorm, "bray") ~ Plant, data = as(sample_data(phyITSnorm), "data.frame"))

#-Example for "sampling period" factor, first growing season
phyITSnorm1 <- subset_samples(phyITSnorm, Growing.Season=="2015-2016") #subset
set.seed(1)
adonis(distance(phyITSnorm1, "bray") ~ Month, data = as(sample_data(phyITSnorm1), "data.frame"))


##----Decomposition of beta-diversity #Require betapart package
#to decompose beta diversity into nestedness / turnover, we use betapart package. 
#Analyses require 2 table of haplotypes presence/absence, composed of the same names for each time comparisons. 
#They are given in the Data-ITS environment

test.t <- beta.temp(Temp1, Temp2, index.family="sor")
test.t
#OR= OilseedRape ; GS1 = First Growing Season; Oct.Dec = comparison between October and December
#beta.sim = Turnover ; beta.sne = Nestedness; beta.sor = total dissimilarity


########Taxonomic composition analyses
##-----Filter genus
#example : Wheat

phyITSnormWheat <- subset_samples(phyITSnorm, Plant=="Wheat")
phyITSbyGenusWheat= tax_glom(phyITSnormWheat, "Genus") #aggregate Genus


#Filter on occurrence
phyITSbyGenusWheatOcc<-phyITSbyGenusWheat
otu_table(phyITSbyGenusWheatOcc)[otu_table(phyITSbyGenusWheatOcc)>=1]<-1 
ConditionTable<-merge_samples(phyITSbyGenusWheatOcc,"Condition")

FilterTable = filter_taxa(ConditionTable, function(x) max(x) > 2, TRUE)
FilterTable
FilterTable_Genus<-get_taxa_unique(FilterTable, taxonomic.rank = "Genus")
keep_taxa<-taxa_names(FilterTable)

GenusTable_Wheat<-prune_taxa(keep_taxa, phyITSbyGenusWheat)
GenusTable_Wheat

#Find and delete Unclassified Genus
Unclassified<-grep("unclassified", FilterTable_Genus, value = TRUE)

GenusTable_Wheat_final <- subset_taxa(GenusTable_Wheat, !Genus %in% Unclassified)
GenusTable_Wheat_final <- subset_taxa(GenusTable_Wheat_final, !Genus %in% "Unclassified")

#Make heatmap on Genera
sample_data(GenusTable_Wheat_final)$Condition<- factor(sample_data(GenusTable_Wheat_final)$Condition, 
                                                       levels = c("Wheat (monoculture).2015-2016.October", "Wheat (monoculture).2015-2016.December", "Wheat (monoculture).2015-2016.February", "Wheat (monoculture).2015-2016.May",
                                                                  "Wheat (monoculture).2016-2017.July", "Wheat (monoculture).2016-2017.October", "Wheat (monoculture).2016-2017.December", "Wheat (monoculture).2016-2017.February", "Wheat (monoculture).2016-2017.May", "Wheat (monoculture).2017-2018.July",
                                                                  "Wheat (Rotation).2015-2016.October", "Wheat (Rotation).2015-2016.December", "Wheat (Rotation).2015-2016.February", "Wheat (Rotation).2015-2016.May",
                                                                  "Wheat (Rotation).2016-2017.July", "Wheat (Rotation).2016-2017.October", "Wheat (Rotation).2016-2017.December", "Wheat (Rotation).2016-2017.February", "Wheat (Rotation).2016-2017.May", "Wheat (Rotation).2017-2018.July"))



plot_heatmap(GenusTable_Wheat_final,taxa.label = "Genus")+facet_grid(~Condition, scales = "free")

