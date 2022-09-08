rm(list = ls())
library(ggplot2)
library(tidyverse)
library(rstatix)
library(qiime2R)

#Location - /home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/

#16s- Taxa bar plot
metadata<-read_q2metadata("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/metadata.NOD_16s.tsv")
SVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/table_16s_filtered_2.qza")$data
taxonomy<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/taxonomy_16s_silva.qza")$data %>% parse_taxonomy()
taxonomy_gg<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/taxonomy_16s_gg.qza")$data %>% parse_taxonomy()


taxasums_genus<-summarize_taxa(SVs, taxonomy)$Genus
taxasums_family<-summarize_taxa(SVs, taxonomy)$Family
taxasums_species<-summarize_taxa(SVs, taxonomy)$Species

taxasums_species_gg<-summarize_taxa(SVs, taxonomy_gg)$Species
taxasums_genus_gg<-summarize_taxa(SVs, taxonomy_gg)$Genus

taxa_barplot(taxasums_genus, metadata, "Type",ntoplot = 10) + ggtitle("Bacterial diversity - Genus level")+theme(plot.title = element_text(hjust = 0.5,size = 11,face="bold"))+
  theme(axis.text=element_text(size=9,colour = "black"), axis.title=element_text(size=11,face="bold"))+ theme(strip.text.x = element_text(size = 9,face="bold")) + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +  theme(legend.text=element_text(size=9)) +
  theme(panel.border = element_rect(size = 2))


ggsave("barplot_genus_16s.pdf",height=12, width=14, device="pdf") # save a PDF 4 inches by 8 inches




taxa_barplot(taxasums_species, metadata, "Type",ntoplot = 10) + ggtitle("Bacterial diversity - Species level")+theme(plot.title = element_text(hjust = 0.5,size = 11,face="bold"))+
  theme(axis.text=element_text(size=9,colour = "black"), axis.title=element_text(size=11,face="bold"))+ theme(strip.text.x = element_text(size = 9,face="bold")) + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +  theme(legend.text=element_text(size=9)) +
  theme(panel.border = element_rect(size = 2))


ggsave("barplot_species_16s.pdf",height=12, width=14, device="pdf") # save a PDF 4 inches by 8 inches


taxa_barplot(taxasums_family, metadata, "Type",ntoplot = 10) + ggtitle("Bacterial diversity - Family level")+theme(plot.title = element_text(hjust = 0.5,size = 9,face="bold"))+
  theme(axis.text=element_text(size=9,colour = "black"), axis.title=element_text(size=11,face="bold"))+ theme(strip.text.x = element_text(size = 9,face="bold")) + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) +  theme(legend.text=element_text(size=9)) +
  theme(panel.border = element_rect(size = 2))



ggsave("barplot_family_16s.pdf", height=12, width=14, device="pdf") # save a PDF 4 inches by 8 inches

##########################################################
#PCoA plots

unweighted_unifrac <-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results-filtered-4/unweighted_unifrac_pcoa_results.qza")
weighted_unifrac <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results-filtered-4/weighted_unifrac_pcoa_results.qza")
bray_curtis <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results-filtered-4/bray_curtis_pcoa_results.qza")



unweighted_unifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Type`, size=10)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + guides(size="none")+
  scale_color_discrete(name="Sample")+
  ggtitle("Bacterial diversity- unweighted unifrac")+theme(plot.title = element_text(hjust = 0.5,size = 20,face="bold"))+  theme(text = element_text(size = 20))

ggsave("PCoA_unweighted_16s.pdf", height=4, width=8, device="pdf") # save a PDF 3 inches by 4 inches

weighted_unifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Type`, size=10)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + guides(size="none")+
  scale_color_discrete(name="Sample")+
  ggtitle("Bacterial diversity- weighted unifrac")+theme(plot.title = element_text(hjust = 0.5,size = 20,face="bold")) +  theme(text = element_text(size = 20)) 




ggsave("PCoA_weighted_16s.pdf", height=4, width=8, device="pdf") # save a PDF 3 inches by 4 inches


bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Type`, size=10)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() + guides(size="none")+
  scale_color_discrete(name="Sample")+
  ggtitle("Bacterial diversity- Bray curtis")+theme(plot.title = element_text(hjust = 0.5,size = 20,face="bold")) +  theme(text = element_text(size = 20)) 


ggsave("PCoA_bray_curtis_16s.pdf", height=4, width=8, device="pdf") # save a PDF 3 inches by 4 inches

###########################################################################################################################
#Alpha diversity statistical tests - Shannon

shannon_filtered_4<-read.csv2("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/shannon_16_filtered_4.tsv",sep = "")


shannon_filtered_4 <- shannon_filtered_4[-1,]

class(shannon_filtered_4$shannon) = "Numeric"
storage.mode(shannon_filtered_4$shannon) = "double"

res.kruskal_shannon_filtered <- kruskal_test(shannon_filtered_4,shannon ~ Type)


summary(res.kruskal_shannon_filtered)  
res.kruskal_shannon_filtered

dunn_test(shannon_filtered_4, shannon ~ Type,p.adjust.method = "bonferroni")


ggplot(shannon_filtered_4,mapping = aes(x=Type,y=shannon))+geom_boxplot()+geom_point()+labs(title ="Shannon diversity - Bacteria (16s)", x = "", y = "Shannon value")+scale_y_continuous(limits=c(0, 3))+theme_q2r() + theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5)) 


mean(filter(shannon_filtered_4,Type == "Kombucha-Day30")[,4])
mean(filter(shannon_filtered_4,Type == "Kombucha-Day1")[,4])
mean(filter(shannon_filtered_4,Type == "Tea-Day30")[,4])
mean(filter(shannon_filtered_4,Type == "Tea-Day1")[,4])


kruskal.test(kombucha_only$shannon ~ kombucha_only$Type)
##############
#Alpha diversity statistical tests - Observed features

ASVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results-filtered-4/observed_features_vector.qza")$data %>% rownames_to_column("SampleID") 

ASVs_1 <- merge(ASVs,metadata,by = "SampleID")

class(ASVs_1$observed_features) = "Numeric"
storage.mode(ASVs_1$observed_features) = "double"

res.kruskal_ASVs <- kruskal_test(ASVs_1,observed_features ~ Type)


summary(res.kruskal_ASVs)  
res.kruskal_ASVs


ggplot(ASVs_1,mapping = aes(x=Type,y=observed_features))+geom_boxplot()+geom_point()+labs(title ="Shannon diversity - Bacteria (16s)", x = "", y = "Observed features")+scale_y_continuous(limits=c(0, 50))+theme_q2r() + theme(text = element_text(size = 20)) + theme(plot.title = element_text(hjust = 0.5)) 

#########################################################################################################################################################################################
