#this script can be use on 16S data 
#Data used for these analyses are given in the R environment 16S_Data.RData

####################Environnement
# Haplo : Haplotypes list with 300 samples. (150 biological samples* 2 technical replicas)
# newname : file used to rename sample names and delete haplotypes not present in the two technical replicas ("Uncertain presences")
# sample.data16S : description of biological samples
# tax_table16S: taxonomic assignement of haplotypes performed on RDP trainset 14 database. In this script, we include the taxonomy table made with : taxa16S <- assignTaxonomy(haplo.table16Sf, "path_to_database", multithread=TRUE)
# Temp1 & Temp2 : haplotypes presence tables. Files used for decomposition of beta-diversity with betapart Package

######## Packages required + custom functions
Packages <- c("phyloseq", "data.table", "ggplot2", "plyr","dplyr","reshape2","grid",
              "gridExtra","scales","dplyr", "ggpubr","vegan","multcompView","rcompanion","betapart")

lapply(Packages, library, character.only = TRUE)


####### "Uncertain presence" removing

Haplo<-Haplo[order(match(rownames(newnames), Haplo[1])), ]
Haplo[1]<-newnames[1]

Haplo[Haplo==0]<-NA
NewHaplo.Table<-aggregate(. ~ X,Haplo, sum, na.action = na.pass)
NewHaplo.Table[is.na(NewHaplo.Table)] <- 0

NewHaplo.Table.with.rownames <- data.frame(NewHaplo.Table[,-1], row.names=NewHaplo.Table[,1])
NewHaplo.Table2<-as.matrix(NewHaplo.Table.with.rownames,row.names=TRUE, header=T)

haplo.table16S<-t(NewHaplo.Table2)
haplo.table16Sf = haplo.table16S[ rowSums(haplo.table16S)!=0, ] 


##############
######## Phyloseq object preparation [phy16S]
tax.table16S <- replace(taxa16S, is.na(taxa16S), "Unclassified")
sample.data16S$Month <- factor(sample.data16S$Month, levels = c("Jul. ", "Oct. ", "Dec. ", "Feb. ","May"))

phy16S <- phyloseq(otu_table(haplo.table16Sf, taxa_are_rows = TRUE),
                 tax_table(tax.table16S),
                 sample_data(sample.data16S))

phy16S

#Removing haplo=0
condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phy16S, condition, 1)
phy16S<-prune_taxa(taxaToKeep, phy16S)
phy16S

#Removing unclassified haplotypes / Anthophyta, Arthropoda, Cercozoa,unclassified_Plantae
table(tax_table(phy16S)[, "Phylum"])
phy16S <- subset_taxa(phy16S, !Phylum %in% c("", "Cyanobacteria/Chloroplast", "Unclassified"))


######## Reads number and normalization
TotalReads = sample_sums(phy16S)

#sequencing depth
SeqDepth = data.table(as(sample_data(phy16S), "data.frame"),
                  TotalReads = sample_sums(phy16S), keep.rownames = TRUE)
pSeqDepth = ggplot(SeqDepth, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
print(pSeqDepth)

#Normalization
phy16Snorm <- transform_sample_counts(phy16S, function(x) round(x/sum(x) *1000000 ))
#write.csv(cbind((otu_table(phy16Snorm)),tax_table(phy16Snorm)), "16S-data.csv")


######## Alpha-diversity analyses : Shannon index
cols <- c("Oilseed Rape" = "#ea0b16", "Wheat (monoculture)" = "#29ce37", "Wheat (Rotation)" = "#007bcf")
measures=c("Shannon")
	#various possible indexes : "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"

RichPlot16S<-plot_richness(phy16Snorm, x="Month", measures = measures)+ 
	stat_compare_means(label.y=c(0.2),label.x.npc = 0.2, method="kruskal.test")+
	geom_boxplot(outlier.shape =NA, aes(fill=Plant.Succession), lwd=0.15, width = 0.3) +#boxplot
	theme_gray(base_size = 10) + ylab(paste(measures,"index (Fungal diversity)")) + xlab("")

RichPlot16S+facet_grid(Growing.Season~Plant.Succession)


##----Wilcoxon rank sum test on Shannon index -- Example for 1st growing season
richness16S<-estimate_richness(phy16Snorm,measures = measures)

test<-sample.data16S[match(rownames(richness16S), rownames(sample.data16S)), ]
compareRich16S<-cbind(test,richness16S)


Year1<-subset(compareRich16S, Growing.Season=="2015-2016")

Oilseed16S1 <- subset(Year1, Plant.Succession=="Oilseed Rape") 
WheatM16S1 <- subset(Year1, Plant.Succession=="Wheat (monoculture)")
WheatR16S1 <- subset(Year1, Plant.Succession=="Wheat (Rotation)")

#example : wilcoxon test between plant succession
compare_means(data=Year1,Shannon~Plant.Succession,method = "wilcox.test")

#example : pairwise wilcoxon for wheat monoculture, 1st growing season
WheatM1.16S.pairwise<-pairwise.wilcox.test(WheatM16S1$Shannon,WheatM16S1$Month, p.adjust.method = "none")
WheatM1.16S.pvalue = fullPTable(WheatM1.16S.pairwise$p.value)
multcompLetters(WheatM1.16S.pvalue) #require : rcompanion package

######## Beta-diversity analyses : Multidimensional scaling on Bray-curtis index
OTU.ord16S<-ordinate(phy16Snorm,"MDS","bray")
Mds16S=plot_ordination(phy16Snorm,OTU.ord16S,type="sample",color="Plant.Succession",shape="Growing.Season")+
	scale_shape_manual(values = c(23,16,17)) + geom_point(size=3)+
	facet_grid(.~Month)+ 
	scale_color_manual(values=cols)+
	theme(panel.spacing = unit(1, "lines"))+
	theme_bw()
Mds16S
##----Permanova test on Bray-curtis index -- Example for "Plant" factor on overall dataset
adonis(distance(phy16Snorm, "bray") ~ Plant, data = as(sample_data(phy16Snorm), "data.frame"))

#-Example for "sampling period" factor, first growing season
phy16Snorm1 <- subset_samples(phy16Snorm, Growing.Season=="2015-2016") #subset
set.seed(1)
adonis(distance(phy16Snorm1, "bray") ~ Month, data = as(sample_data(phy16Snorm1), "data.frame"))


##----Decomposition of beta-diversity #Require betapart package
#to decompose beta diversity into nestedness / turnover, we use betapart package. 
#Analyses require 2 table of haplotypes presence/absence, composed of the same names for each time comparisons. 
#They are given in the environnment XXX

beta.div.dissimilarity <- beta.temp(Temp1, Temp2, index.family="sor")
beta.div.dissimilarity
#OR= OilseedRape ; GS1 = First Growing Season; Oct.Dec = comparison between October and December
#beta.sim = Turnover ; beta.sne = Nestedness; beta.sor = total dissimilarity


########Taxonomic composition analyses
##-----Filter genus
#example : Wheat

phy16SnormWheat <- subset_samples(phy16Snorm, Plant=="Wheat")
phy16SbyGenusWheat= tax_glom(phy16SnormWheat, "Genus") #aggregate Genus


#Filter on occurrence
phy16SbyGenusWheatOcc<-phy16SbyGenusWheat
otu_table(phy16SbyGenusWheatOcc)[otu_table(phy16SbyGenusWheatOcc)>=1]<-1 
ConditionTable<-merge_samples(phy16SbyGenusWheatOcc,"Condition")

FilterTable = filter_taxa(ConditionTable, function(x) max(x) > 2, TRUE)
FilterTable
FilterTable_Genus<-get_taxa_unique(FilterTable, taxonomic.rank = "Genus")
keep_taxa<-taxa_names(FilterTable)

GenusTable_Wheat<-prune_taxa(keep_taxa, phy16SbyGenusWheat)
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

