# how I made the toy dataset, in case anybody (future me) needs to know

require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('~/Projects/rcc/sceMNN')

# get TCR contigs
sample_ids <- system('ls ~/Projects/rcc/Data/TCR/samples', intern = TRUE)
contigs <- NULL
for (ii in 1:length(sample_ids)) {
    tcr_ii <- read.csv(paste0("~/Projects/rcc/Data/TCR/samples/", sample_ids[ii],
                              "/filtered_contig_annotations.csv"))
    tcr_ii$barcode <- gsub("1$", ii, tcr_ii$barcode)
    tcr_ii$type <- "TCR"
    tcr_ii$clonotype_tcr <- paste0(sample_ids[ii], "-", tcr_ii$raw_clonotype_id)
    tcr_ii$sample <- sample_ids[ii]
    contigs <- rbind(contigs, tcr_ii)
}
rm(ii, tcr_ii)

tcr.list <- split(DataFrame(contigs), factor(contigs$barcode, sce$cell))
class(tcr.list)

sce$contigs <- tcr.list

# let's just look at T cells
toy <- sce[,sce$clusTCELL != -1]

# and narrow it down to one patient
toy <- toy[,toy$sample %in% c('S11_N','S11_T')]

# and downsample a bit
idx <- c(sample(which(toy$sample=='S11_N'), 500),
         sample(which(toy$sample=='S11_T'), 500))
toy <- toy[,idx]

# # convert counts to sparse matrix, so no messing with HDF5
require(Matrix)
assay(toy, 'counts') <- as(assay(toy, 'counts'), 'dgCMatrix')

# make new SCE with just the parts we want
toy2 <- SingleCellExperiment(assays = SimpleList(counts = assay(toy,'counts')), 
                             colData = colData(toy), 
                             reducedDims = list(TSNE = reducedDim(toy,'TSNE'), PCAMNN = reducedDim(toy,'corrected')))
toy <- toy2

saveRDS(toy, file='~/Desktop/toyTCRdata.rds')


