
# add TCR (or similar) data to SCE

dir <- '~/Projects/OLD/rcc/Data/TCR/samples'
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('~/Projects/OLD/rcc/sceMNN')

# future work: make it more flexible so that 'sce' can be a SCE object or a directory,
#    use DropletUtils::read10xCounts to get SCE from directory

# need to figure out barcodes. Why are they all *_1 in the TCR data? is replacing that with ii really doing the matching properly? (94.4%)

require(S4Vectors)
addTCRtoSCE <- function(dir, sce){
    # get TCR contigs
    sample_ids <- system(paste('ls',dir), intern = TRUE)
    contigs <- NULL
    for (ii in 1:length(sample_ids)) {
        tcr_ii <- read.csv(file.path(dir, sample_ids[ii],
                                  "filtered_contig_annotations.csv"))
        tcr_ii$barcode <- gsub("1$", ii, tcr_ii$barcode)
        tcr_ii$type <- "TCR"
        tcr_ii$clonotype_tcr <- paste0(sample_ids[ii], "-", tcr_ii$raw_clonotype_id)
        tcr_ii$sample <- sample_ids[ii]
        contigs <- rbind(contigs, tcr_ii)
    }
    rm(ii, tcr_ii)
    # convert to 
    tcr.list <- split(DataFrame(contigs), factor(contigs$barcode, sce$cell))
    class(tcr.list)
    
    sce$contigs <- tcr.list
    return(sce)
}


#sce <- addTCRtoSCE(dir, sce)
# inspection (look for mismatched samples)
# tcrsamp <- sapply(tcr.list[,'sample'], function(x){ ifelse(length(x)>0, x[1], NA)})
# tcrsamp[is.na(tcrsamp)] <- sce$sample[is.na(tcrsamp)]
# x <- sce[, tcrsamp != sce$sample]
