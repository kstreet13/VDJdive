# source('EMquant.R')
# sce <- readRDS('~/Desktop/toyTCRdata.rds')
# sce <- EMquant(sce, sample = 'sample')


# cell-level counts by sample
#' @param sce SingleCellExperiment object with clonotype counts (sce$clono)
#' @param by variable (name or vector) by which to split clonotype counts

splitClonotypes <- function(sce, by){
    if(is.null(sce$clono)){
        stop('No clonotype counts found.')
    }
    if(length(by) == 1){
        byVar <- factor(sce[[by]])
    }else{
        stopifnot(length(by) == ncol(sce))
        byVar <- factor(by)
    }
    out <- lapply(sort(unique(byVar)), function(lv){
        sce$clono[which(byVar == lv), ,drop = FALSE]
    })
    names(out) <- as.character(sort(unique(byVar)))
    return(out)
}

# get sample level counts
#' @param sce SingleCellExperiment object with clonotype counts (sce$clono)
#' @param by variable (name or vector) by which to summarize clonotype counts

summarizeClonotypes <- function(sce, by){
    if(is.null(sce$clono)){
        stop('No clonotype counts found.')
    }
    if(length(by) == 1){
        byVar <- factor(sce[[by]])
    }else{
        stopifnot(length(by) == ncol(sce))
        byVar <- factor(by)
    }
    out <- t(sapply(sort(unique(byVar)), function(lv){
        colSums(sce$clono[which(byVar == lv), ,drop = FALSE])
    }, USE.NAMES = TRUE))
    rownames(out) <- as.character(sort(unique(byVar)))
    return(out)
}







