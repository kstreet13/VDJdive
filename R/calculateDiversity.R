
# diversity functions
.totalCount <- function(p){
    return(sum(p))
}
.shannon <- function(p){
    p <- p[p > 0]
    p <- p / sum(p)
    return(-sum(p * log(p)))
}
.normentropy <- function(p){
    p <- p[p > 0]
    p <- p / sum(p)
    return(-sum(p * log(p)) / log(length(p)))
}
.invsimpson <- function(p){
    p <- p / sum(p)
    return(1 / sum(p^2))
}
.ginisimpson <- function(p){
    p <- p / sum(p)
    return(1 - sum(p^2))
}
.chao1 <- function(f){
    # f is the table of singletons, doubletons, etc.
    # for integer counts, if p = x/sum(x), then f = table(x)
    return(sum(f) + f[1]*(f[1]-1) / (2*(f[2]+1)))
}
#' @importFrom stats qnorm
.chaobunge <- function(f, t = 10, conf = 0.95){
    # adapted from jipingw/SPECIES::ChaoBunge()
    
    # f is the table of singletons, doubletons, etc.
    # for integer counts, if p = x/sum(x), then f = table(x)
    f <- f[seq_len(max(which(f>0)))] # trim trailing 0s
    if (t!=round(t)||t<0){
        stop("The cutoff t to define less abundant species must be", 
             " non-negative integer!")
    } 
    if(is.numeric(conf)==FALSE||conf>1||conf<0){
        stop("confidence level must be a numerical value between 0 and 1,",
             " e.g. 0.95")
    } 
    if(t > length(f)){
        warning("The t that defines the abundant/rare species must be no ", 
                "larger than the most abundant species!","\n",
                "We use t=", length(f), " in this calculation!")
        t <- length(f)
    }
    m <- length(f)

    ############################
    A1 <- sum(f[seq_len(t)])
    A2 <- sum(f[seq_len(t)])-f[1]
    B1 <- sum(seq_len(t) * f[seq_len(t)])
    #B2 <- B1 - f[1]
    #C2 <- sum(seq_len(t) * (seq_len(t)-1) * f[seq_len(t)])
    
    ## point estimate (Equation 2, page 533 of Chao and Bunge 2002)
    if(t == length(f)){
        theta <- 1 - f[1]*sum(seq_len(m)^2*f) / (sum(seq_len(m)*f))^2
        chao40 <- sum(f[seq_len(t)[-1]]) / theta
        chao4 <- (1/theta-1) * sum(f[seq_len(m)[-1]]) - f[1]+sum(f)
    }
    if(t < length(f)){
        theta <- 1 - f[1]*sum(seq_len(t)^2*f[seq_len(t)]) /
            (sum(seq_len(t)*f[seq_len(t)]))^2
        chao4 <- sum(f[t+seq_len(length(f)-t)]) +
            sum(f[seq_len(t)[-1]]/theta)
        chao40 <- sum(f[seq_len(t)[-1]])/theta
    }
    
    ## duplication standard error(Unnumbered Eq. right below Eq. 2, 
    ## page 534 of Chao and Bunge 2002)
    D <- sum(seq_len(t)*seq_len(t)*f[seq_len(t)])
    partial3 <- numeric(t)
    partial3[1] <- sum(-f[seq_len(t)[-1]] / (1-f[1]*D/B1^2)^2 *
                           (-(D+f[1])/(B1^2) + (2*f[1]*D/B1^3)))
    
    for (i in seq_len(t)[-1]){
        partial3[i] <- 1/(1-f[1]*D/B1^2) - A2/(1-f[1]*D/B1^2)^2 *
            (-f[1]*i^2/B1^2+2*f[1]*D*i/B1^3)
    }
    
    cova3 <- matrix(0, nrow = t, ncol = t)
    
    for (i in seq_len(t)){
        for (j in seq_len(t)){
            if(i==j){
                cova3[i,j] <- f[i]*(1-f[i]/chao40)		
            }
            if(i!=j){
                cova3[i,j] <- -f[i]*f[j]/chao40	
            }	
        }
    }
    
    SE <- 0
    for (i in seq_len(t)){
        for (j in seq_len(t)){
            SE <- SE + partial3[i] * partial3[j] * cova3[i,j]
        }   
    }
    SE <- partial3 %*% cova3 %*% partial3
    SE4 <- sqrt(SE[1,1])
    
    ##confidence interval
    coe <- qnorm(1-(1-conf)/2, 0, 1)
    C <- exp(coe*log(1+SE4^2/(chao4-A1)^2)^0.5)
    lb <- floor(A1+(chao4-A1)/C)
    ub <- ceiling(A1+(chao4-A1)*C)
    # CI0 <- matrix(c(lb,ub),1,2)
    # colnames(CI0) <- c("lb","ub")
    # return(list(est=chao4, CI=CI0))
    return(c(est = chao4, CI.lower = lb, CI.upper = ub))
}


#' @title Sample diversity estimation
#' @name calculateDiversity
#' @param ... Additional arguments passed to external calculation methods.
#' @import methods
#' @export
setGeneric(name = "calculateDiversity",
           signature = "x",
           def = function(x, ...) standardGeneric("calculateDiversity"))


#' @rdname calculateDiversity
#'
#' @description This function uses various methods to estimate the clonotypic
#'   diversity of samples based on a matrix of clonotype abundances (samples are
#'   columns).
#'
#' @param x A matrix of abundance values where rows are features (clonotypes)
#'   and columns are samples. This is created with \code{summarizeClonotypes}
#'   using a sparse matrix computed with either \code{EMquant} or
#'   \code{CRquant}.
#' @param methods A character vector specifying which diversity measures to use
#'   (default = \code{'all'}, see Details).
#'
#' @details Available methods are total count (\code{'count'}), Shannon entropy
#'   (\code{'shannon'}), Simpson index (\code{'simpson'}), inverse Simpson index
#'   (\code{'invsimpson'}), Chao1 richness (\code{'chao1'}), and Chao-Bunge
#'   richness (\code{'chaobunge'}). A special value of \code{'all'} is also
#'   allowed, which will run all methods listed above.
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
#' data('contigs')
#' x <- clonoStats(contigs)
#' calculateDiversity(x)
#'
#' @export
setMethod(f = "calculateDiversity",
          signature = signature(x = "list"),
          definition = function(x, methods = c('all','count','shannon',
                                               'normentropy', 'invsimpson',
                                               'ginisimpson', 'chao1',
                                               'chaobunge'),
                                ...){
              methods <- match.arg(methods, several.ok = TRUE)
              if('all' %in% methods){
                  methods <- c('count', 'shannon','normentropy','invsimpson',
                               'ginisimpson','chao1','chaobunge')
              }

              # check if all counts are integers
              ints <- all(x$abundance %% 1 == 0)
              if(!ints & any(c('chao1','chaobunge') %in% methods)){
                  warning('Methods "chao1" and "chaobunge" are not valid with ',
                          'non-integer abundances.')
              }
              
              p <- x$abundance
              f <- x$frequency
              
              # loop over methods
              results <- lapply(methods, function(m){
                  if(m == 'count'){
                      return(vapply(seq_len(ncol(p)), function(j){
                          .totalCount(p[,j])
                      }, FUN.VALUE = 0.0))
                  }
                  if(m == 'shannon'){
                      return(vapply(seq_len(ncol(p)), function(j){
                          .shannon(p[,j])
                      }, FUN.VALUE = 0.0))
                  }
                  if(m == 'normentropy'){
                      return(vapply(seq_len(ncol(p)), function(j){
                          .normentropy(p[,j])
                      }, FUN.VALUE = 0.0))
                  }
                  if(m == 'invsimpson'){
                      return(vapply(seq_len(ncol(p)), function(j){
                          .invsimpson(p[,j])
                      }, FUN.VALUE = 0.0))
                  }
                  if(m == 'ginisimpson'){
                      return(vapply(seq_len(ncol(p)), function(j){
                          .ginisimpson(p[,j])
                      }, FUN.VALUE = 0.0))
                  }
                  if(m == 'chao1'){
                      return(vapply(seq_len(ncol(f)), function(j){
                          .chao1(f[-1,j]) # remove the 0s row
                      }, FUN.VALUE = 0.0))
                  }
                  if(m == 'chaobunge'){
                      cb <- vapply(seq_len(ncol(f)), function(j){
                          .chaobunge(f[-1,j]) # remove the 0s row
                      }, FUN.VALUE = rep(0.0,3))
                      rownames(cb) <- c('chaobunge.est','chaobunge.CI.lower',
                                        'chaobunge.CI.upper')
                      return(cb)
                  }
              })
              names(results) <- methods
              results <- do.call(rbind, results)
              colnames(results) <- colnames(p)
              
              return(results)
          })


