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
#' @param x Matrix or SingleCellExperiment with clonotype counts (x$clono)
#' @param by variable (name or vector) by which to summarize clonotype counts
#' @param clonoCol when x is an SCE, the name for the clonotype matrix in the colData (default = 'clono')
#' @return A table of sample-level clonotype counts
require(Matrix)

setGeneric(name = "summarizeClonotypes",
           signature = c("x","by"),
           def = function(x, by, ...) standardGeneric("summarizeClonotypes"))

setMethod(f = "summarizeClonotypes",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, by, clonoCol = 'clono'){
              if(is.null(x[[clonoCol]])){
                  stop('No clonotype counts found.')
              }
              if(length(by) == 1){
                  byVar <- factor(x[[by]])
              }else{
                  byVar <- factor(by)
              }
              out <- summarizeClonotypes(x[[clonoCol]], byVar)
              return(out)
          })

setMethod(f = "summarizeClonotypes",
          signature = signature(x = "matrix"),
          definition = function(x, by){
              summarizeClonotypes(Matrix(x), by)
          })

setMethod(f = "summarizeClonotypes",
          signature = signature(x = "Matrix"),
          definition = function(x, by){
              stopifnot(length(by) == nrow(x))
              by <- factor(by)
              out <- t(sapply(sort(unique(by)), function(lv){
                  colSums(x[which(by == lv), ,drop = FALSE])
              }, USE.NAMES = TRUE))
              rownames(out) <- as.character(sort(unique(by)))
              return(out)
          })





