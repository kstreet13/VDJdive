
# cell-level counts by sample
#' @title Split cell-level clonotype counts by sample
#' @name splitClonotypes
#' @param ... additional arguments.
#' @export
setGeneric(name = "splitClonotypes",
           signature = c("x","by"),
           def = function(x, by, ...) standardGeneric("splitClonotypes"))

#' @rdname splitClonotypes
#'
#' @description Takes a matrix of cell-level clonotype counts and splits them
#'   into a list of group-specific counts (typically samples).
#'
#' @param x A \code{Matrix} of cell-level clonotype assignments
#'   (cells-by-clonotypes) or a \code{SingleCellExperiment} object with such a
#'   matrix stored in the \code{clono} slot of the \code{colData}.
#' @param by A character vector or factor by which to split the clonotype
#'   counts. If \code{x} is a \code{SingleCellExperiment} object, this can also
#'   be a character, giving the name of the column from the \code{colData} to
#'   use as this variable.
#'
#' @return A list of \code{Matrix} objects providing the cell-level counts for
#'   each unique value of \code{by} (if \code{by} denotes sample labels, each
#'   matrix in the list will contain the cells from a single sample).
#'
#' @examples
#' example(addVDJtoSCE)
#' x <- clonoStats(contigs, assignment = TRUE)
#' splitClonotypes(x$assignment, by = sce$sample)
#'
#' @importClassesFrom Matrix Matrix
#' @export
setMethod(f = "splitClonotypes",
          signature = signature(x = "Matrix"),
          definition = function(x, by){
              stopifnot(length(by) == nrow(x))
              if(!is.factor(by)){
                  by <- factor(by)
              }
              out <- lapply(levels(by), function(lv){
                  x[which(by == lv), ,drop = FALSE]
              })
              names(out) <- as.character(levels(by))
              return(out)
          })


#' @rdname splitClonotypes
#' @importFrom Matrix Matrix
#' @export
setMethod(f = "splitClonotypes",
          signature = signature(x = "matrix"),
          definition = function(x, by){
              splitClonotypes(Matrix(x), by)
          })

#' @rdname splitClonotypes
#' @param clonoStats character. If \code{x} is a \code{SingleCellExperiment},
#'   the name of the element in the \code{metadata} of \code{x} that contains
#'   the output of \code{clonoStats}. Must include cell-level clonotype
#'   assignments (ie. \code{assignment = TRUE}).
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
setMethod(f = "splitClonotypes",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, by, clonoStats = 'clonoStats'){
              if(is.null(metadata(x)[[clonoStats]])){
                  stop('No clonotype counts found.')
              }
              if(length(by) == 1){
                  byVar <- factor(x[[by]])
              }else{
                  byVar <- factor(by)
              }
              return(splitClonotypes(metadata(x)[[clonoStats]]$assignment, 
                                     byVar))
          })





#' @title Get sample-level clonotype counts
#' @name summarizeClonotypes
#' @param ... additional arguments.
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
#' @param clonoStats character. If \code{x} is a \code{SingleCellExperiment},
#'   the name of the element in the \code{metadata} of \code{x} that contains
#'   the output of \code{clonoStats}. Must include cell-level clonotype
#'   assignments (ie. \code{assignment = TRUE}).
#' @param mode Type of summarization to perform. Default is \code{'sum'}, which
#'   sums clonotype abundances within each sample (or level of \code{'by'}).
#'   Alternative is \code{'tab'}, which constructs a table of clonotype
#'   frequencies (ie. singletons, doubletons, etc.) by sample.
#'
#' @return A matrix clonotype counts where each row corresponds to a unique
#'   value of \code{by} (if \code{by} denotes sample labels, this is a matrix of
#'   sample-level clonotype counts).
#'
#' @examples
#' example(addVDJtoSCE)
#' x <- clonoStats(contigs, assignment = TRUE)
#' summarizeClonotypes(x$assignment, by = sce$sample)
#'
#' @importClassesFrom Matrix Matrix
#' @importFrom Matrix rowSums
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "Matrix"),
          definition = function(x, by, mode = c('sum','tab')){
              mode <- match.arg(mode)
              if(!is.factor(by)){
                  by <- factor(by)
              }
              stopifnot(length(by) == nrow(x))
              if(mode == 'sum'){
                  out <- vapply(levels(by), function(lv){
                      colSums(x[which(by == lv), ,drop = FALSE])
                  }, FUN.VALUE = rep(0,ncol(x)))
              }
              if(mode == 'tab'){
                  lim <- max(table(by))+1
                  if(all(x@x %% 1 == 0)){ # integer counts
                      out <- vapply(levels(by), function(lv){
                          ind <- which(by==lv)
                          return(table(factor(colSums(x[ind, ,drop=FALSE]),
                                              levels = seq_len(lim)-1)))
                      }, FUN.VALUE = rep(0,lim))
                  }else{ # non-integer counts
                      out <- vapply(levels(by), function(lv){
                          sub.x <- x[which(by==lv), ,drop=FALSE]
                          # I tried putting the following code in a function,
                          # but it didn't work for some mysterious reason:
                          # Error in intI(i, n = x@Dim[1], dn[[1]], give.dn = FALSE) : 
                          #     index larger than maximal 1
                          # from "p.dists[2, ind1] <- EX[ind1]"
                          nz <- unname(Matrix::colSums(sub.x > 0)) # non-zero count
                          EX <- unname(Matrix::colSums(sub.x)) # sum of probs
                          p.dists <- Matrix(0, nrow = lim, ncol = ncol(sub.x), 
                                            sparse = TRUE)
                          # absent clonotypes
                          ind0 <- unname(which(nz == 0))
                          p.dists[1, ind0] <- 1
                          # clonotypes that could only be in one cell
                          ind1 <- unname(which(nz == 1))
                          p.dists[1, ind1] <- 1 - EX[ind1]
                          p.dists[2, ind1] <- EX[ind1]
                          # clonotypees that could only be in (up to) two cells
                          ind2 <- unname(which(nz == 2))
                          P2 <- matrix(sub.x[, ind2, drop = FALSE]@x, nrow = 2)
                          p.dists[1, ind2] <- exp(colSums(log(1-P2)))
                          p.dists[3, ind2] <- exp(colSums(log(P2)))
                          p.dists[2, ind2] <- 1 - p.dists[1, ind2] - 
                              p.dists[3, ind2]
                          # then do the rest the slow way
                          for(jj in which(nz > 2)){
                              probs <- sub.x[, jj, drop=FALSE]@x
                              distr <- 1
                              for(p in probs){
                                  distr <- c(distr*(1-p), 0) + c(0, distr*p)
                              }
                              p.dists[,jj] <- c(distr, 
                                                rep(0, lim-length(distr)))
                          }
                          return(Matrix::rowSums(p.dists))
                      }, FUN.VALUE = rep(0,lim))
                  }
                  # trim excess 0s
                  out <- out[seq_len(max(c(0,which(rowSums(out) > 0)))), ,
                             drop = FALSE]
                  rownames(out) <- seq_len(nrow(out))-1
              }
              return(out)
          })

#' @rdname summarizeClonotypes
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, by = "sample", 
                                clonoStats = 'clonoStats', ...){
              if(is.null(metadata(x)[[clonoStats]])){
                  stop('No clonotype counts found.')
              }
              if(length(by) == 1){
                  byVar <- factor(x[[by]])
              }else{
                  byVar <- factor(by)
              }
              return(summarizeClonotypes(metadata(x)[[clonoStats]]$assignment,
                                         byVar, ...))
          })

#' @rdname summarizeClonotypes
#' @importFrom Matrix Matrix
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "matrix"),
          definition = function(x, by, ...){
              summarizeClonotypes(Matrix(x), by, ...)
          })

