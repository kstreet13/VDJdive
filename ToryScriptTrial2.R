#ToyScriptTrial2
library(SingleCellExperiment)
library(dplyr)
toyTCRdata <- readRDS("/Users/mercedeh/Desktop/PostDoColabProj/toyTCRdata.rds")
#MakeFile
# obtain coldata matrix as data.frame (rownames are group_name)
tcr_coldata <- data.frame(attributes(toyTCRdata)[["colData"]][, 1:21])
tcr_coldata$group_name <- rownames(tcr_coldata)
# obtain contigs matrix as data.frame, joined with tcr_coldata
tcr_matrix <- inner_join(data.frame(toyTCRdata$contigs), tcr_coldata,
                         by=c("group_name", "sample"),
                         suffix=c("_contigs", "_coldata"))#deal with duplicate columns
# filter
tcr_file<-dplyr::filter(tcr_matrix, cdr3 != "None" & productive == "True") #Filtering cdr3 is None values and full length false values
tcr_file_fin<-dplyr::filter(tcr_file, full_length == "True")


