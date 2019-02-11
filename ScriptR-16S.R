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

sum(otu_table(phy16S))

#Removing haplo=0
condition <- function(x) { x > 0 } 
taxaToKeep <- genefilter_sample(phy16S, condition, 1)
phy16S<-prune_taxa(taxaToKeep, phy16S)
phy16S

#Removing unclassified haplotypes / Anthophyta, Arthropoda, Cercozoa,unclassified_Plantae
table(tax_table(phy16S)[, "Phylum"])
phy16S <- subset_taxa(phy16S, !Phylum %in% c("", "Cyanobacteria/Chloroplast", "Unclassified"))


phy16S<-merge_phyloseq(phy16S,tree)
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
##

######## Alpha-diversity analyses : Shannon index
cols <- c("Oilseed Rape" = "#ea0b16", "Wheat (monoculture)" = "#29ce37", "Wheat (Rotation)" = "#007bcf")
measures=c("Shannon")
	#various possible indexes : "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"

RichPlot16S<-plot_richness(phy16Snorm, x="Month", measures = measures)+ 
	stat_compare_means(label.y=c(0.2),label.x.npc = 0.2, method="kruskal.test")+
	geom_boxplot(outlier.shape =NA, aes(fill=Plant.Succession), lwd=0.15, width = 0.3) +#boxplot
	theme_gray(base_size = 10) + ylab(paste(measures,"index (Fungal diversity)")) + xlab("")

RichPlot16S+facet_grid(Growing.Season~Plant.Succession)



##----Alpha div
richness16S<-estimate_richness(phy16Snorm,measures = measures)

test<-sample.data16S[match(rownames(richness16S), rownames(sample.data16S)), ]
compareRich16S<-cbind(test,richness16S)

#Comparison Wheat in monoculture - wheat in rotation
Wheat<-subset(compareRich16S, Plant.Succession != "Oilseed Rape")
model<-lm(Shannon~Plant.Succession + Growing.Season+ Month.ent  , data=Wheat)
library(car)
Anova(model)


#Comparison Oilseed Rape- Wheat in rotation
Rotation<-subset(compareRich16S, Plant.Succession != "Wheat (monoculture)")
model<-lm(Shannon~Plant.Succession + Growing.Season+ Month.ent+Parcelle.Code  , data=Rotation)
Anova(model)

####### Faith's Diveristy ######

PhyPD<-subset_samples(phy16S, Saison !="S3")
library(picante)
PhyPDmatrix<-t(as(otu_table(PhyPD),"matrix"))
FD<-pd(samp = PhyPDmatrix, tree = phy_tree(PhyPD), include.root = F)
FD<-cbind(FD,sample_data(PhyPD))


RichFD16S<- ggplot(FD, aes(x=Month, y=PD))+
  geom_boxplot(aes(fill=Plant.Succession))+stat_compare_means(label.y =40,method="kruskal.test")+
  facet_grid(Growing.Season~Plant.Succession, scales="free_x")

RichFD16S+theme_bw()+ylim(0,42)

Year1<-subset(FD, Growing.Season=="2016-2017")


######## Beta-diversity analyses : Multidimensional scaling on Bray-curtis index
OTU.ord16S<-ordinate(phy16Snorm,"MDS","bray")
Mds16S=plot_ordination(phy16Snorm,OTU.ord16S,type="sample",color="Plant.Succession",shape="Growing.Season")+
	scale_shape_manual(values = c(23,16,17)) + geom_point(size=3)+
	facet_grid(.~Month)+ 
	scale_color_manual(values=cols)+
	theme(panel.spacing = unit(1, "lines"))+
	theme_bw()
Mds16S


OTU.ord16S<-ordinate(PhywithTree,"MDS","UniFrac")
Mds16S=plot_ordination(phy16Snorm,OTU.ord16S,type="sample",color="Plant.Succession",shape="Growing.Season")+
  scale_shape_manual(values = c(23,16,17)) + geom_point(size=3)+
  facet_grid(.~Month)+ 
  scale_color_manual(values=cols)+
  theme(panel.spacing = unit(1, "lines"))+
  theme_bw()
Mds16S

##----Permanova test on Bray-curtis index -- Resoumission
Wheatbray<-subset_samples(phy16Snorm,Plant.Succession !="Oilseed Rape")
set.seed(1)
adonis2(distance(Wheatbray, "bray") ~ Plant.Succession+Growing.Season+Month.ent,by="margin", data = as(sample_data(Wheatbray), "data.frame"))

Rotationbray<-subset_samples(phy16Snorm,Plant.Succession !="Wheat (monoculture)")
set.seed(1)
adonis2(distance(Rotationbray, "bray") ~ Plant.Succession+Growing.Season+Month.ent+Parcelle.Code,by="margin", data = as(sample_data(Rotationbray), "data.frame"))



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

phy16SbyGenus= tax_glom(phy16Snorm, "Genus") #aggregate Genus


#Filter on occurrence
phy16SbyGenusOcc<-phy16SbyGenus
otu_table(phy16SbyGenusOcc)[otu_table(phy16SbyGenusOcc)>=1]<-1 
ConditionTable<-merge_samples(phy16SbyGenusOcc,"Condition")

FilterTable = filter_taxa(ConditionTable, function(x) max(x) > 2, TRUE)
FilterTable
FilterTable_Genus<-get_taxa_unique(FilterTable, taxonomic.rank = "Genus")
keep_taxa<-taxa_names(FilterTable)

GenusTable_Wheat<-prune_taxa(keep_taxa, phy16SbyGenus)
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

plot_heatmap(ent10,  method="MDS",  taxa.order = taxa_names(nameX),taxa.label = "Genus", na.value = "#FFFFFF")+
  facet_grid(~Plant.Succession+Growing.Season+Month, scales="free")+
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        strip.background = element_rect(color = "black", size = 0.5))+
  theme(axis.text.y = element_text(face="italic",size=10))+
  scale_fill_gradient2(low="#ffffeb", mid=	"green", high="#000000",  
                       na.value = "white", trans = log_trans(4), 
                       midpoint = log(16000, base = 4))


