library(data.table)
library(BiocManager)
library(decontam)

# Import features
Bacterial <- read.csv("feature-table_16s.tsv",sep ="")
Fungal <- read.csv("feature-table_ITS.tsv",sep ="")

# Filter control samples
Controls_bacteria <- subset(Bacterial,select = c("Tea.control.1","Tea.control.2","Water.control.1","Water.control.2"))
Controls_bacteria_0 <- subset(Controls_bacteria,Tea.control.1 == 0 || Tea.control.2 == 0 || Water.control.1 == 0 ||  Water.control.2 == 0)


Bacterial <- as.data.frame(t(as.matrix(Bacterial)))
Fungal <- as.data.frame(t(as.matrix(Fungal)))


Bacterial_data <- as.matrix((Bacterial))
Fungal_data <- as.matrix((Fungal))

# Include a vector indicate the negative controls
neg <- c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)

# Find contaminant bacteria 
Contaminant_bacteria <- isContaminant(Bacterial_data,neg = neg,detailed = TRUE,threshold = 0.5)
Contaminant_fungi <- isContaminant(Fungal_data,neg = neg,detailed = TRUE,threshold = 0.5)


write.csv(Contaminant_bacteria,file = "Contaminants_16s",sep = " ")
write.csv(Contaminant_fungi,file = "Contaminants_ITS",sep = " ")

Contaminant_bacteria <- read.csv("Contaminants_16s",sep = ",")
Contaminant_fungi <- read.csv("Contaminants_ITS",sep = ",")

Filterd_bac <- subset(Contaminant_bacteria,contaminant == "TRUE")
Filterd_fungi <- subset(Contaminant_fungi,contaminant == "TRUE")

Filterd_bac <- subset(Filterd_bac,select=-c(freq,prev,p.freq,p.prev,p,contaminant))
Filterd_fungi <- subset(Filterd_fungi,select=-c(freq,prev,p.freq,p.prev,p,contaminant))

names(Filterd_bac)[names(Filterd_bac) == "X"] <- "feature-id"
names(Filterd_fungi)[names(Filterd_fungi) == "X"] <- "feature-id"

# Save the contaminant features to use it in qiime to remove from the actual data before analysis
write.table(Filterd_bac,file = "Contam_16s.txt",sep = " ",row.names = FALSE,quote = FALSE)
write.table(Filterd_fungi,file = "Contam_ITS.txt",sep = " ",row.names = FALSE,quote = FALSE)

