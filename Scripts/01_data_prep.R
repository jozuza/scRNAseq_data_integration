##R script created on November, 2021
# This script creates an indivudual seurat object for each dataset 
# Work from your 'Script' dir

#Load libraries
library(readxl)
library(useful)
library(dplyr)
library(Seurat)

# Gene annotation
Genes.Ann <- read.table("../external_metadata/GeneId_Symbol.txt", header = T, sep = "\t")
Genes.Ann$ID_Symbol = paste(Genes.Ann$EnsemblID, Genes.Ann$Symbol, sep = "_")

#-------------------------------------------------------------------------------
# Xiang et al, 2020 data--------------------------------------------------------
# Load GSE136447 metadadata
# From supplementary table8
metasheet <- read.table("../external_metadata/GSE136447/Embryo_Info.txt",
                                header = T, sep = "\t")
names(metasheet) = c(
  "Sample", "Day", "Embryo", "lineage"
)
metasheet$Cell = 1:nrow(metasheet)

count <- read.table("../RawData/GSE136447/GSE136447_555sample_gene_count_matrix.txt.gz",
                    header = T, sep = "\t")

# genes at the PAR region are duplicated in the count matrix. There are 45 of them
count = count[grep("PAR_Y",count$gene_id, invert = T),]
count$gene_id <-  sapply(strsplit(count$gene_id, "\\."), `[`, 1)
rownames(count) = count$gene_id
count$gene_id <- NULL


count$EnsemblID = rownames(count)
count = merge(count, Genes.Ann, by = "EnsemblID")
count = unique(count)
rownames(count) = count$ID_Symbol
count$EnsemblID <- NULL
count$Symbol <- NULL
count$ID_Symbol <- NULL

# Get the same cells (Samples) from metasheet and count matrix
metasheet = metasheet %>% filter(Sample %in% colnames(count)); dim(metasheet)
# 555    5
count = count[,which(colnames(count)%in% metasheet$Sample)]; dim(count)
# 64792   555
#put them in the same order
count <- count[,match(metasheet$Sample, colnames(count))]
metasheet <- metasheet[order(match(metasheet$Sample,colnames(count))),]

#Check point
identical(colnames(count), metasheet$Sample)
# [1] TRUE

#create_Seurat ----------------------------------------------------
library(Seurat)
DataSet <- CreateSeuratObject(
  counts = count, #a matrix-like object with unnormalized data with cells as columns and features as rows
  project = "GSE136447", #Project name 
  min.cells = 3, #Include features detected in at least this many cells. Will subset the counts matrix as well.
  min.features = 25 #Include cells where at least this many features are detected.
)

DataSet
DataSet@meta.data$Sample  <- metasheet$Sample
DataSet@meta.data$Day     <- metasheet$Day
DataSet@meta.data$Embryo  <- metasheet$Embryo
DataSet@meta.data$lineage <- metasheet$lineage
DataSet@meta.data$Cell    <- metasheet$Cell

#Check point
identical(rownames(DataSet@meta.data),DataSet@meta.data$Sample)
# [1] TRUE

#save -------------------------------------------------------------
write.table(count, 
            file="../RawData/GSE136447/GSE136447_count.txt", quote=F,col.names=T,row.names=T, sep="\t")

write.table(metasheet, "../external_metadata/GSE136447/GSE136447_metadata.txt", 
            quote=F, col.names=T,row.names=F, sep="\t")

saveRDS(DataSet, file = "../RawData/GSE136447/GSE136447_SeuratObj.rds")

rm(count, Samples, metasheet, DataSet)
#-------------------------------------------------------------------------------
# Zhou et al, 2019 full lenght data---------------------------------------------
# Load HRR000128 metadadata
# 
metasheet <- read.table("../external_metadata/HRR000128/HRR000128_Metadata.txt",
                        header = T, sep = "\t")
metasheet$Sample.Name <- NULL
metasheet$Embryo.Name <- NULL

names(metasheet) = c(
  "Sample", "lineage", "Sex","Day", "Embryo"
)
metasheet$Cell = 1:nrow(metasheet)

count <- read.table("../RawData/HRA000128/HRA000128_read_count.txt",
                    header = T, sep = "\t")

#Add gene names
count$EnsemblID = count$Gene_ID
count$Gene_ID <- NULL
count = merge(count, Genes.Ann, by = "EnsemblID")
count = unique(count)
rownames(count) = count$ID_Symbol
count$EnsemblID <- NULL
count$ID_Symbol <- NULL

# Get the same cells (Samples) from metasheet and count matrix
metasheet = metasheet %>% filter(Sample %in% colnames(count)); dim(metasheet)
# 233   6
count = count[,which(colnames(count)%in% metasheet$Sample)]; dim(count)
# 60675   233
#put them in teh same order
count <- count[,match(metasheet$Sample, colnames(count))]
metasheet <- metasheet[order(match(metasheet$Sample,colnames(count))),]

#Check point
identical(colnames(count), metasheet$Sample)
# [1] TRUE

#create_Seurat ----------------------------------------------------
library(Seurat)
DataSet <- CreateSeuratObject(
  counts = count, #a matrix-like object with unnormalized data with cells as columns and features as rows
  project = "HRR000128", #Project name 
  min.cells = 3, #Include features detected in at least this many cells. Will subset the counts matrix as well.
  min.features = 25 #Include cells where at least this many features are detected.
)

DataSet
DataSet@meta.data$Sample  <- metasheet$Sample
DataSet@meta.data$Day     <- metasheet$Day
DataSet@meta.data$Embryo  <- metasheet$Embryo
DataSet@meta.data$Sex     <- metasheet$Sex
DataSet@meta.data$lineage <- metasheet$lineage
DataSet@meta.data$Cell    <- metasheet$Cell

#Check point
identical(rownames(DataSet@meta.data),DataSet@meta.data$Sample)
# [1] TRUE

#save -------------------------------------------------------------
saveRDS(DataSet, file = "../RawData/HRA000128/HRA000128_SeuratObj.rds")

rm(count, metasheet, DataSet)


# Petropoulos et al, 2013 data--------------------------------------------------
# Load EMTAB3929 data
emtab3929_meta <- readRDS("../RawData/EMTAB3929/EMTAB3929.rds")
# Get count matrix
count = emtab3929_meta@ExperimentList@listData$gene@assays$data$count

# (1) Remove version numbers from Ensembl gene IDs

rownames(count) <- sapply(strsplit(rownames(count), "\\."), `[`, 1)
dim(count) # 65218  1529

# check embryo stages and remove E3 and E4
Samples = colnames(count)
unique(gsub(Samples, pattern="?[.].*", replacement=""))
# "E3" "E4" "E5" "E6" "E7"

Samples = Samples[grep("E3", Samples, invert = T)]
Samples = Samples[grep("E4", Samples, invert = T)]
unique(gsub(Samples, pattern="?[.].*", replacement=""))

# keep only embryos with the 3 lineages represented 

# (3) Cell lineage information. Compile sample metasheet.
cell_lineage_data <- read_xlsx("../external_metadata/EMTAB3929/stirparo2018_tableS4.xlsx", sheet = 1)

cell_lineage_data <- cell_lineage_data[cell_lineage_data$Study == "Petropoulos et al., 2016 (ERP012552)", ] # 1,481 cells
cell_lineage_data$Cell <- gsub("_", ".", cell_lineage_data$Cell)
cell_lineage_data$EStage <- sapply(strsplit(cell_lineage_data$Embryo, "_"), "[", 1)

metasheet = cell_lineage_data %>% select("Embryo","Cell","Stage","EStage",
                                         "Revised lineage (this study)")

names(metasheet) =c("Embryo","Sample","Stage", "Day", "lineage")
metasheet$Cell <- sapply(strsplit(metasheet$Sample, "\\."), tail, n = 1)
metasheet = metasheet %>% filter(Sample %in% Samples)
rm(Samples)

count$EnsemblID <- NULL
count$Symbol <- NULL
count <- na.omit(count)
rownames(count) = count$ID_Symbol

#Add gene names
count$EnsemblID = rownames(count)
count = merge(count, Genes.Ann, by = "EnsemblID")
count = unique(count)
rownames(count) = count$ID_Symbol


# Get the same cells (Samples) from metasheet and count matrix
metasheet = metasheet %>% filter(Sample %in% colnames(count)); dim(metasheet)
# 1218    9
count = count[,which(colnames(count)%in% metasheet$Sample)]; dim(count)
# 65218   1218

#put them in teh same order
count <- count[,match(metasheet$Sample, colnames(count))]
metasheet <- metasheet[order(match(metasheet$Sample,colnames(count))),]

#check point
identical(colnames(count), metasheet$Sample)
# [1] TRUE

# Save objects for integration and analysis

write.table(count, 
            file="../RawData/EMTAB3929/EMTAB3929_count.txt", quote=F,col.names=T,row.names=T, sep="\t")

write.table(metasheet, "../external_metadata/EMTAB3929/EMTAB3929_metadata.txt", 
            quote=F, col.names=T,row.names=F, sep="\t")

#create_Seurat ----------------------------------------------------
DataSet <- CreateSeuratObject(
  counts = count, #a matrix-like object with unnormalized data with cells as columns and features as rows
  project = "EMTAB3929", #Project name 
  min.cells = 3, #Include features detected in at least this many cells. Will subset the counts matrix as well.
  min.features = 25 #Include cells where at least this many features are detected.
)

DataSet
DataSet@meta.data$Sample  <- metasheet$Sample
DataSet@meta.data$Day     <- metasheet$Day
DataSet@meta.data$Embryo  <- metasheet$Embryo
DataSet@meta.data$lineage <- metasheet$lineage
DataSet@meta.data$Cell    <- metasheet$Cell
DataSet@meta.data$Stage    <- metasheet$Stage

#check point
identical(rownames(DataSet@meta.data),DataSet@meta.data$Sample)
# [1] TRUE

#save -------------------------------------------------------------
saveRDS(DataSet, file = "../RawData/EMTAB3929/EMTAB3929_SeuratObj.rds")

rm(emtab3929_meta, count, Samples, cell_lineage_data, metasheet, DataSet)

#-------------------------------------------------------------------------------


