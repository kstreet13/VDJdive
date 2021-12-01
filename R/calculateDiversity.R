# devtools::install_github('adw96/breakaway')

# usage:
# sce <- readRDS('~/Desktop/toyTCRdata.rds')
# sce <- EMquant(sce, sample = 'sample')
# k <- t(summarizeClonotypes(sce, 'sample'))
# calculateDiversity(k)


.div_function <- function(k, method, ints, scale_factor, ...) {
    if (method %in% c("simpson", "shannon", "invsimpson")) {
        # simpson, shannon and invsimpson calculated using diversity() from vegan package
        ret <- data.frame(t(vegan::diversity(t(k), ..., index=method)))
    } else if (method == "chao1") {
        # use fossil::chao1() for chao1
        if(!ints){
            k <- round(k * scale_factor)
        }
        ret <- data.frame(t(apply(k, 2, fossil::chao1, ...))) #apply chao1 function for every column
        if(!ints){
            ret <- ret / scale_factor
        }
    } else if (method == "chaobunge") {
        # use breakaway() for chaobunge
        if(!ints){
            k <- round(k * scale_factor)
        }
        ret <- data.frame(apply(k, 2, function(x, ...) {
            bw <- breakaway::breakaway(x, ...)
            c(estimate=bw$estimate, error=bw$error)
        }))
        if(!ints){
            ret <- ret / scale_factor
        }
    } 
    return(ret)
}


#' @title Sample diversity estimation
#' @name calculateDiversity
#' @param ... Additional arguments passed to external calculation methods.
#' @export
setGeneric(name = "calculateDiversity",
           signature = "k",
           def = function(k, ...) standardGeneric("calculateDiversity"))


#' @rdname calculateDiversity
#' 
#' @description This function uses various methods to estimate the clonotypic
#'   diversity of samples based on a matrix of clonotype abundances (samples are
#'   columns).
#' 
#' @param k A matrix of abundance values where rows are features (clonotypes)
#'   and columns are samples.
#' @param methods A character vector specifying which diversity measures to use
#'   (default = \code{'all'}, see Details).
#' @param scale_factor Numeric scaling factor for handling non-integer counts
#'   (see Details).
#' 
#' @details Available methods are Shannon entropy (\code{'shannon'}), Simpson
#'   index (\code{'simpson'}), inverse Simpson index (\code{'invsimpson'}),
#'   Chao1 richness (\code{'chao1'}), and Chao-Bunge richness
#'   (\code{'chaobunge'}). A special value of \code{'all'} is also allowed,
#'   which will run all methods listed above.
#' 
#' @details The \code{'chao1'} and \code{'chaobunge'} estimates assume all
#'   abundances are integers. When this is not the case for the input matrix,
#'   \code{k}, all values are multiplied by the \code{scaling_factor} and
#'   rounded to the nearest integer. The resulting estimate is then divided by
#'   \code{scaling_factor} to return to the original scale. The
#'   \code{'shannon'}, \code{'simpson'}, and \code{'invsimpson'} methods work
#'   with any input type.
#'   
#' @return A matrix of diversity estimates for each sample. Note that the
#'   \code{'chaobunge'} method also includes an estimate of the standard error.
#'
#' @examples 
#' data('example_contigs')
#' samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
#' counts <- EMquant(contigs)
#' k <- t(summarizeClonotypes(counts, samples))
#' calculateDiversity(k)
#'
#' @importFrom vegan diversity
#' @importFrom fossil chao1
#' @importFrom breakaway breakaway
#' 
#' @export
setMethod(f = "calculateDiversity",
          signature = signature(k = "matrix"),
          definition = function(k, methods = c('all','shannon','simpson',
                                               'invsimpson','chao1',
                                               'chaobunge'), 
                                scale_factor = 100, ...){
              methods <- match.arg(methods, several.ok = TRUE)
              if('all' %in% methods){
                  methods <- c('shannon','simpson','invsimpson','chao1','chaobunge')
              }
              # any parm but k or methods is ... for e.g. cutoff for breakaway which will pass down
              # function to calculate diversity for a given method
              
              # check if all counts are integers
              ints <- all(k %% 1 == 0)
              
              # loop over methods
              results <- sapply(methods, function(method){
                  .div_function(k, method, ints, scale_factor, ...)
              })
              colnames(results) <- paste0(colnames(results), '.estimate')
              if('chaobunge' %in% methods){
                  ests <- sapply(results[,'chaobunge.estimate'], function(x){x[1]})
                  errs <- sapply(results[,'chaobunge.estimate'], function(x){x[2]})
                  cbres <- cbind(chaobunge.estimate = ests, chaobunge.error = errs)
                  results <- results[,colnames(results) != 'chaobunge.estimate']
                  results <- cbind(results, cbres)
              }
              
              return(t(results))
          })


