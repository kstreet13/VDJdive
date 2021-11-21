
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




#' @title Get sample-level clonotype counts
#' @name summarizeClonotypes
#' @export
setGeneric(name = "summarizeClonotypes",
           signature = c("x","by"),
           def = function(x, by, ...) standardGeneric("summarizeClonotypes"))

#' @rdname summarizeClonotypes
#' 
#' @description Takes a matrix of cell-level clonotype counts and sums them
#'   within groups (typically samples).
#' @param x A (usually sparse) matrix of cell-level clonotype counts (cells are
#'   rows and clonotypes are columns). Alternatively, a
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} with such a matrix
#'   stored in the \code{colData}.
#' @param by A character vector or factor by which to summarize the clonotype
#'   counts. If \code{x} is a \code{SingleCellExperiment} object, this can also
#'   be a character, giving the name of the column from the \code{colData} to
#'   use as this variable.
#' @param clonoCol A character providing the name of the cell-level clonotype
#'   counts matrix in the \code{colData} of \code{x} (default = \code{'clono'}).
#'   Only applies when \code{x} is a \code{SingleCellExperiment} object.
#'   
#' @return A matrix clonotype counts where each row corresponds to a unique
#'   value of \code{by} (if \code{by} denotes sample labels, this is a matrix of
#'   sample-level clonotype counts).
#' 
#' @export
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

#' @rdname summarizeClonotypes
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
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

#' @rdname summarizeClonotypes
#' @import Matrix
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "matrix"),
          definition = function(x, by){
              summarizeClonotypes(Matrix(x), by)
          })






