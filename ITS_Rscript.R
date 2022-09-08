rm(list = ls())
library(ggplot2)
library(tidyverse)
library(qiime2R)

#Location - /home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/

#ITS - taxa bar plot
metadata<-read_q2metadata("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/metadata.NOD_ITS.tsv")
SVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/table_ITS_filtered_2.qza")$data
taxonomy<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/taxonomy_ITS.qza")$data %>% parse_taxonomy()

taxasums_genus<-summarize_taxa(SVs, taxonomy)$Genus
taxasums_family<-summarize_taxa(SVs, taxonomy)$Family
taxasums_species<-summarize_taxa(SVs, taxonomy)$Species

taxa_barplot(taxasums_genus, metadata, "Type",ntoplot = 10) + ggtitle("Fungal diversity - Genus level")+theme(plot.title = element_text(hjust = 0.5,size = 11,face="bold"))+
  theme(axis.text=element_text(size=9,colour = "black"), axis.title=element_text(size=11,face="bold"))+ theme(strip.text.x = element_text(size = 9,face="bold")) + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +  theme(legend.text=element_text(size=9)) +
  theme(panel.border = element_rect(size = 2))


ggsave("barplot_genus_ITS.pdf", height=12, width=14, device="pdf") 

taxa_barplot(taxasums_species, metadata, "Type",ntoplot = 10) + ggtitle("Fungal diversity - Species level")+theme(plot.title = element_text(hjust = 0.5,size = 11,face="bold"))+
  theme(axis.text=element_text(size=9,colour = "black"), axis.title=element_text(size=11,face="bold"))+ theme(strip.text.x = element_text(size = 9,face="bold")) + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +  theme(legend.text=element_text(size=9)) +
  theme(panel.border = element_rect(size = 2))

ggsave("barplot_species_ITS.pdf", height=12, width=14, device="pdf") # save a PDF 4 inches by 8 inches


taxa_barplot(taxasums_family, metadata, "Type",ntoplot = 10) + ggtitle("Fungal diversity - Family level")+theme(plot.title = element_text(hjust = 0.5,size = 9,face="bold"))+
  theme(axis.text=element_text(size=9,colour = "black"), axis.title=element_text(size=11,face="bold"))+ theme(strip.text.x = element_text(size = 9,face="bold")) + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +  theme(legend.text=element_text(size=9)) +
  theme(panel.border = element_rect(size = 2))



ggsave("barplot_family_ITS.pdf", height=12, width=14, device="pdf") # save a PDF 4 inches by 8 inches

###########################################################################################
#PcoA plots

unweighted_unifrac <-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/unweighted_unifrac_pcoa_results.qza")
weighted_unifrac <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/weighted_unifrac_pcoa_results.qza")
bray_curtis <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/bray_curtis_pcoa_results.qza")


unweighted_unifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Type`, size=10)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + guides(size="none")+
  scale_color_discrete(name="Sample")+
  ggtitle("Fungal diversity- unweighted unifrac")+theme(plot.title = element_text(hjust = 0.5,size = 20,face="bold")) +  theme(text = element_text(size = 20)) 

ggsave("PCoA_unweighted_ITS.pdf", height=4, width=8, device="pdf") # save a PDF 3 inches by 4 inches



weighted_unifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Type`,size=10)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + guides(size="none")+
  scale_color_discrete(name="Sample")+
  ggtitle("Fungal diversity- weighted unifrac")+theme(plot.title = element_text(hjust = 0.5,size = 20,face="bold")) +  theme(text = element_text(size = 20)) 


ggsave("PCoA_weighted_ITS.pdf", height=4, width=8, device="pdf") # save a PDF 3 inches by 4 inches


bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Type`, size=10)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + guides(size="none")+
  scale_color_discrete(name="Sample")+
  ggtitle("Fungal diversity- Bray curtis")+theme(plot.title = element_text(hjust = 0.5,size = 20,face="bold")) +  theme(text = element_text(size = 20)) 


ggsave("PCoA_bray_curtis_ITS.pdf", height=4, width=8, device="pdf") # save a PDF 3 inches by 4 inches



################################################################################
#Alpha diversity statistical tests - Shannon

shannon_filtered_2<-read.csv2("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/shannon_ITS_filtered_2.tsv",sep = "")


shannon_filtered_2 <- shannon_filtered_2[-1,]

class(shannon_filtered_2$shannon) = "Numeric"
storage.mode(shannon_filtered_2$shannon) = "double"

res.kruskal_shannon_filtered_2 <- kruskal_test(shannon_filtered_2,shannon ~ Type)


res.kruskal_shannon_filtered_2

dunn_test(shannon_filtered_2, shannon ~ Type,p.adjust.method = "bonferroni")

ggplot(shannon_filtered_2,mapping = aes(x=Type,y=shannon))+geom_boxplot()+geom_point()+labs(title =" Shannon diversity - Fungal (ITS) ", x = "", y = "Shannon value")+scale_y_continuous(limits=c(0, 3))+theme_q2r() + theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5))


#########
#Alpha diversity statistical tests - Observed features


ASVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/observed_otus_vector.qza")$data %>% rownames_to_column("SampleID") 

ASVs_1 <- merge(ASVs,metadata,by = "SampleID")

class(ASVs_1$observed_otus) = "Numeric"
storage.mode(ASVs_1$observed_otus) = "double"

res.kruskal_ASVs <- kruskal_test(ASVs_1,observed_otus ~ Type)


summary(res.kruskal_ASVs)  
res.kruskal_ASVs

dunn_test(ASVs_1,observed_otus ~ Type,p.adjust.method = "bonferroni")


ggplot(ASVs_1,mapping = aes(x=Type,y=observed_otus))+geom_boxplot()+geom_point()+labs(title ="Shannon diversity - Bacteria (16s)", x = "", y = "Observed features")+scale_y_continuous(limits=c(0, 50))+theme_q2r() + theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5)) 
########################################################################################
