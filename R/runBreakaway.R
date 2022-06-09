#' @title Clonotype richness estimation with Breakaway
#' @name runBreakaway
#'
#' @description This function uses the Breakaway method to estimate the
#'   clonotype richness (total number of clonotypes) present in each group of a
#'   \code{clonoStats} object.
#'
#' @param x A \code{\link{clonoStats}} object.
#' @param nof1 logical. Indicates whether to use the \code{breakaway_nof1}
#'   function, for abundance data that may contain spurious singletons.
#' @param ... Additional arguments passed to \code{breakaway} or
#'   \code{breakaway_nof1}.
#'
#' @return A list of \code{alpha_estimate} objects, one per group, containing
#'   detailed results of running the Breakaway estimator on the vector of
#'   clonotype frequencies from that group.
#'
#' @references 
#' Willis, A. and Bunge, J. (2015). Estimating diversity via frequency ratios.
#' \emph{Biometrics}.
#'
#' Willis, A. (2015). Species richness estimation with high diversity but
#' spurious singletons. \emph{arXiv}.
#'
#' @examples
#' data('contigs')
#' x <- clonoStats(contigs, method = 'unique')
#' runBreakaway(x)
#'
#' @export
runBreakaway <- function(x, nof1 = FALSE, ...){
    if(requireNamespace('breakaway', quietly = TRUE)){
        nc <- ncol(clonoFrequency(x))
        groupnames <- colnames(clonoFrequency(x))
        if(!nof1){
            out <- lapply(seq_len(nc), function(i){
                breakaway::breakaway(clonoFrequency(x)[-1,i], ...)
            })
        }else{
            out <- lapply(seq_len(nc), function(i){
                breakaway::breakaway_nof1(clonoFrequency(x)[-1,i], ...)
            })
        }
        names(out) <- groupnames
        return(out)
    }else{
        stop('Package "breakaway" must be installed.')
    }
}
