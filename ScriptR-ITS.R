#this script can be use on ITS data 
#Data used for these analyses are given in the R environment ITS_Data.RData

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
sum(otu_table(phyITS))
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


##----Wilcoxon rank sum test on Shannon index 
richnessITS<-estimate_richness(phyITSnorm,measures = measures)

test<-sample.dataITS[match(rownames(richnessITS), rownames(sample.dataITS)), ]
compareRichITS<-cbind(test,richnessITS)

#comparison wheat in monoculture - wheat in rotation
Wheat<-subset(compareRichITS, Plant.Succession != "Oilseed Rape")
model<-lm(Shannon~Plant.Succession + Growing.Season+ Month.ent  , data=Wheat)
library(car)
Anova(model)


#Comparison Wheat in rotation - oilseed rape
Rotation<-subset(compareRichITS, Plant.Succession != "Wheat (monoculture)")
model<-lm(Shannon~Plant.Succession + Growing.Season+ Month.ent+Parcelle.Code  , data=Rotation)
library(car)
Anova(model)


######## Beta-diversity analyses : Multidimensional scaling on Bray-curtis index
OTU.ordITS<-ordinate(phyITSnorm,"MDS","bray")
MdsITS=plot_ordination(phyITSnorm,OTU.ordITS,type="sample",color="Plant.Succession",shape="Growing.Season")+
	scale_shape_manual(values = c(23,16,17)) + geom_point(size=3)+
	facet_grid(.~Month)+ 
	scale_color_manual(values=cols)+
	theme(panel.spacing = unit(1, "lines"))+
	theme_bw()
MdsITS

##----Permanova test on Bray-curtis index 
Wheatbray<-subset_samples(phyITSnorm,Plant.Succession !="Oilseed Rape")
set.seed(1)
adonis2(distance(Wheatbray, "bray") ~ Plant.Succession+Growing.Season+Month.ent,by="margin", data = as(sample_data(Wheatbray), "data.frame"))

Rotationbray<-subset_samples(phyITSnorm,Plant.Succession !="Wheat (monoculture)")
set.seed(1)
adonis2(distance(Rotationbray, "bray") ~ Plant.Succession+Growing.Season+Month.ent+Parcelle.Code,by="margin", data = as(sample_data(Rotationbray), "data.frame"))



##----Decomposition of beta-diversity #Require betapart package
#to decompose beta diversity into nestedness / turnover, we use betapart package. 
#Analyses require 2 table of haplotypes presence/absence, composed of the same names for each time comparisons. 
#They are given in the environnment XXX

test.t <- beta.temp(Temp1, Temp2, index.family="sor")
test.t
#OR= OilseedRape ; GS1 = First Growing Season; Oct.Dec = comparison between October and December
#beta.sim = Turnover ; beta.sne = Nestedness; beta.sor = total dissimilarity


########Taxonomic composition analyses
##-----Filter genus
#example : Wheat

phyITSnormWheat <- phyITSnorm
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
sample_data(GenusTable_Wheat_final)$Plant.Succession<- factor(sample_data(GenusTable_Wheat_final)$Plant.Succession, 
                                                       levels = c("Wheat (monoculture)", "Wheat (Rotation)","Oilseed Rape"))


tax_table(GenusTable_Wheat_final)[,1] <- gsub("k__", "", tax_table(GenusTable_Wheat_final)[,1]);tax_table(GenusTable_Wheat_final)[,2] <- gsub("p__", "", tax_table(GenusTable_Wheat_final)[,2]);tax_table(GenusTable_Wheat_final)[,3] <- gsub("c__", "", tax_table(GenusTable_Wheat_final)[,3]);tax_table(GenusTable_Wheat_final)[,4] <- gsub("o__", "", tax_table(GenusTable_Wheat_final)[,4]);tax_table(GenusTable_Wheat_final)[,5] <- gsub("f__", "", tax_table(GenusTable_Wheat_final)[,5]);tax_table(GenusTable_Wheat_final)[,6] <- gsub("g__", "", tax_table(GenusTable_Wheat_final)[,6])


TopNOTUs <- names(sort(taxa_sums(GenusTable_Wheat_final), TRUE)[1:60])
ent10   <- prune_taxa(TopNOTUs, GenusTable_Wheat_final)
nameX<- rev(tax_table(ent10)[,"Genus"][order(tax_table(ent10)[,"Genus"]),])

plot_heatmap(ent10,  method="MDS", taxa.order = taxa_names(nameX), taxa.label = "Genus", na.value = "#FFFFFF")+
  facet_grid(~Plant.Succession+Growing.Season+Month, scales="free")+
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        strip.background = element_rect(color = "black", size = 0.5))+
  theme(axis.text.y = element_text(face="italic",size=10))+
  scale_fill_gradient2(low="#ffffeb", mid=	"green", high="#000000",  
                       na.value = "white", trans = log_trans(4), 
                       midpoint = log(16000, base = 4))

