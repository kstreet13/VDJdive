#' @importClassesFrom Matrix dgCMatrix
setClassUnion("dgCMatrixOrNULL",members=c("dgCMatrix", "NULL"))

#' @title \code{clonoStats} object class
#' @aliases clonoStats-class
#' 
#' @description The \code{clonoStats} class is designed to hold the output of
#'   the \code{\link{clonoStats}} function. This always includes two group-level
#'   summaries: clonotype abundances and clonotype frequencies. "Group" most
#'   often refers to sample of origin, but may alternatively refer to any
#'   partitioning of cells, such as clusters. Clonotype names are stored
#'   efficiently in two factors. Additionally, a large, sparse matrix of each
#'   cell's clonotype assignment may be included.
#' 
#' @slot abundance Summary table of clonotype abundances within each group
#'   (sample). Provides specific information about how often each clonotype is
#'   observed in each group.
#' @slot frequency Summary table of clonotype frequencies within each group
#'   (sample). Provides a summary of clonotype abundances (ie. number of
#'   singletons, doubletons, etc.).
#' @slot group Factor variable giving the group of origin for each cell.
#' @slot assignment Optional matrix of clonotype assignments in each individual
#'   cell. Rows correspond to cells and row sums are all 0 or 1. When using
#'   \code{clonoStats} with \code{method = "unique"} or \code{method =
#'   "CellRanger"}, all values will be 0 or 1. With \code{method = "EM"},
#'   fractional values are allowed, representing the assignment confidence.
#' @slot names1 Factor listing the first component of each clonotype name (alpha
#'   chains, for TCR data).
#' @slot names2 Factor listing the second component of each clonotype name (beta
#'   chains, for TCR data).
#'   
#' @return An object of class \code{clonoStats}.
#' 
#' @seealso \code{\link{clonoStats}}
#' 
#' @examples 
#' data('contigs')
#' cs <- clonoStats(contigs)
#' cs
#' 
#' @importClassesFrom Matrix dgCMatrix
#' @import methods
#' @export
setClass(
    Class = "clonoStats",
    slots = list(
        abundance = "dgCMatrix",
        frequency = "dgCMatrix",
        group = 'factor',
        assignment = "dgCMatrixOrNULL",
        names1 = "factor",
        names2 = "factor"
    )
)


setValidity("clonoStats", function(object) {
    if(ncol(object@abundance) != ncol(object@frequency)) {
        return(paste("Abundance and frequency matrices must have same number",
                     "of columns"))
    }
    if(!is.null(rownames(object@abundance))){
        return("Abundance matrix must not have rownames")
    }
    if(!is.null(object@assignment)){
        if(ncol(object@assignment) != nrow(object@abundance)){
            return(paste("Assignment matrix must have as many columns as",
                         "abundance matrix has rows"))
        }
        if(!is.null(colnames(object@assignment))){
            return("Assignment matrix must not have colnames")
        }
        if(nrow(object@assignment) != length(object@group)){
            return("Group variable length must match number of rows of",
                   "assignment matrix")
        }
    }
    if(length(object@names1) != nrow(object@abundance)){
        return(paste("First clonotype name component (names1) must have length",
                     "equal to nrow(abundance)"))
    }
    if(length(object@names2) != nrow(object@abundance)){
        return(paste("Second clonotype name component (names2) must have",
                     "length equal to nrow(abundance)"))
    }
    return(TRUE)
})

