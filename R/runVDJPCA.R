
#' @title Run PCA on clonotype abundance matrix
#' @param ... additional arguments.
#' @name runVDJPCA
#' @export
setGeneric(name = "runVDJPCA",
           signature = "k",
           def = function(k, ...) standardGeneric("runVDJPCA"))


#' @rdname runVDJPCA
#' 
#' @description Perform Principal Components Analysis (PCA) on the matrix of
#'   sample-level clonotype abundances. In the context of clonotype analysis,
#'   this is a form of beta diversity.
#' 
#' @param k A matrix of abundance values where rows are features (clonotypes)
#'   and columns are samples.
#' @param unit Character value indicating whether the unit of interest is
#'   \code{"samples"} or \code{"clonotypes"}.
#'   
#' @return A list with class \code{"prcomp"}. The component \code{x} stores the
#'   reduced-dimensional representation of the data. For a full description, see
#'   \code{\link[stats]{prcomp}}.
#' 
#' @examples 
#' data('contigs')
#' x <- clonoStats(contigs)
#' runVDJPCA(x$abundance)
#' 
#' @import stats
#' @export
setMethod(f = "runVDJPCA",
          signature = signature(k = "matrix"),
          definition = function(k, unit = c("samples","clonotypes")){
    unit <- match.arg(unit)
    
    if(unit == "clonotypes"){
        PCA_ret <- prcomp(t(k))
    } else {
        PCA_ret <- prcomp(k) 
    }
    return(PCA_ret)
})
