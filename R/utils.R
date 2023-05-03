#' @include clonoStats_class.R basilisk.R
NULL

#' @importFrom reticulate import
#' @import Matrix
.nonInt_tab <- function(x, lim, lang = c("r","python")){
    lang <- match.arg(lang)
    nz <- unname(colSums(x > 0)) # non-zero count
    EX <- unname(colSums(x)) # sum of probs
    if(lim >= 1){
        # clonotypes that could only be in zero or one cell
        ind1 <- which(nz == 1)
        # skip straight to the row sums
        RS1 <- c(sum(nz == 0) + length(ind1) - sum(EX[ind1]), 
                 sum(EX[ind1]), rep(0, lim-2))
    }else{
        RS1 <- rep(0,lim)
    }
    if(lim >= 2){
        # clonotypes that could only be in (up to) two cells
        ind2 <- which(nz == 2)
        P2 <- matrix(x[, ind2, drop = FALSE]@x, nrow = 2)
        RS2 <- c(sum(exp(colSums(log(1-P2)))), 0,
                 sum(exp(colSums(log(P2)))), rep(0, lim-3))
        RS2[2] <- length(ind2) - RS2[1] - RS2[3]
    }else{
        RS2 <- rep(0, lim)
    }
    if(lim >= 3){
        # then do the rest (either in python or R)
        ind3p <- which(nz >= 3)
        counts <- x[, ind3p, drop = FALSE]@x
        probs_list <- unname(split(counts, 
                                   factor(rep(seq_len(length(ind3p)), 
                                              times = nz[ind3p]))))
        if(lang == 'python'){
            cl <- basiliskStart(pyenv)
            distrs_list <- basiliskRun(proc = cl, function(probs_list){
                mod <- import(module = "vdjHelpers", 
                              convert = TRUE)
                return(mod$make_distrs(probs_list))
            }, probs_list = probs_list)
            basiliskStop(cl)
        }
        if(lang == 'r'){
            # this might cause memory problems
            distrs_list <- lapply(probs_list, function(probs){
                distr <- 1
                for(p in probs){
                    distr <- c(distr*(1-p), 0) + c(0, distr*p)
                }
                return(distr)
            })
        }
        tmpI <- as.integer(unlist(lapply(lengths(distrs_list), function(l){ seq_len(l)-1 })))
        tmpP <- as.integer(c(0,cumsum(lengths(distrs_list))))
        tmp <- new('dgCMatrix', i = tmpI, p = tmpP, x = as.numeric(unlist(distrs_list)),
                   Dim = as.integer(c(lim, length(distrs_list))))
        RS3 <- rowSums(tmp)
    }else{
        RS3 <- rep(0, lim)
    }
    return((RS1+RS2+RS3)[-1])
}

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
#'   use as this variable. Similar to the \code{group} argument for
#'   \code{\link{clonoStats}}.
#' @param contigs character. If \code{x} is a \code{SingleCellExperiment}, the
#'   name of the \code{SplitDataFrameList} in the \code{colData} of \code{x}
#'   containing contig information.
#' @param clonoStats character. If \code{x} is a \code{SingleCellExperiment},
#'   the name of the element in the \code{metadata} of \code{x} that contains
#'   the output of \code{clonoStats}. Must include cell-level clonotype
#'   assignments (ie. \code{assignment = TRUE}).
#' @param mode Type of summarization to perform. Default is \code{'sum'}, which
#'   sums clonotype abundances within each sample (or level of \code{'by'}).
#'   Alternative is \code{'tab'}, which constructs a table of clonotype
#'   frequencies (ie. singletons, doubletons, etc.) by sample.
#' @param lang Indicates which implementation of the \code{"tab"} summarization
#'   to use. Options are \code{'r'} (default) or \code{'python'}. Only used if
#'   non-integer clonotype abundances are present and \code{mode = "tab"}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the
#'   parallel backend for distributed clonotype assignment operations (split by
#'   \code{group}). Default is \code{BiocParallel::SerialParam()}.
#'
#' @return If \code{mode = 'sum'}, returns a matrix clonotype abundances where
#'   each row corresponds to a clonotype and each column a value of \code{by}
#'   (if \code{by} denotes sample labels, this is a matrix of sample-level
#'   clonotype counts). If \code{mode = 'tab'}, returns a matrix of clonotype
#'   frequencies, where each row corresponds to a frequency (singletons,
#'   doubletons, etc.) and each column a value of \code{by}.
#'
#' @examples
#' example(addVDJtoSCE)
#' x <- clonoStats(contigs, assignment = TRUE)
#' summarizeClonotypes(x, by = sce$sample)
#'
#' @import Matrix
#' @importFrom BiocParallel SerialParam bplapply
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "Matrix"),
          definition = function(x, by, mode = c('sum','tab'), 
                                lang = c('r','python'),
                                BPPARAM = SerialParam()){
              mode <- match.arg(mode)
              lang <- match.arg(lang)
              if(!is.factor(by)){
                  by <- factor(by)
              }
              stopifnot(length(by) == nrow(x))
              if(mode == 'sum'){
                  out <- bplapply(levels(by), function(lv){
                      colSums(x[which(by == lv), ,drop = FALSE])
                  }, BPPARAM = BPPARAM)
                  out <- do.call(cbind, out)
              }
              if(mode == 'tab'){
                  lim <- max(table(by))
                  if(all(x@x %% 1 == 0)){ # integer counts
                      out <- bplapply(levels(by), function(lv){
                          ind <- which(by==lv)
                          return(table(factor(colSums(x[ind, ,drop=FALSE]),
                                              levels = seq_len(lim))))
                      }, BPPARAM = BPPARAM)
                  }else{ # non-integer counts
                      out <- bplapply(levels(by), function(lv){
                          sub.x <- x[which(by==lv), ,drop=FALSE]
                          return(.nonInt_tab(sub.x, lim, lang))
                      }, BPPARAM = BPPARAM)
                  }
                  out <- do.call(cbind, out)
                  
                  # trim excess 0s
                  out <- out[seq_len(max(c(0,which(rowSums(out) > 0)))), ,
                             drop = FALSE]
                  rownames(out) <- seq_len(nrow(out))
              }
              colnames(out) <- levels(by)
              return(Matrix(out, sparse = TRUE, doDiag = FALSE))
          })

#' @rdname summarizeClonotypes
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, by = "sample", contigs = 'contigs',
                                clonoStats = 'clonoStats', ...){
              if(is.null(metadata(x)[[clonoStats]])){
                  stop('No clonotype counts found.')
              }
              if(length(by) == 1){ # check SCE and contigs
                  byVar <- factor(x[[by]]) # SCE
                  if(length(byVar) == 0){
                      if(by %in% names(x[[contigs]][[1]])){
                          byVar <- factor(vapply(seq_len(ncol(x)), function(i){
                              x[[contigs]][,by][[i]][1]
                          }, FUN.VALUE = 'a')) # contigs
                      }
                  }
              }else{
                  byVar <- factor(by)
              }
              return(summarizeClonotypes(clonoAssignment(
                  metadata(x)[[clonoStats]]), byVar, ...))
          })

#' @rdname summarizeClonotypes
#' @import Matrix
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "matrix"),
          definition = function(x, by, ...){
              summarizeClonotypes(Matrix(x), by, ...)
          })

#' @rdname summarizeClonotypes
#' @export
setMethod(f = "summarizeClonotypes",
          signature = signature(x = "clonoStats"),
          definition = function(x, by, ...){
              if(is.null(clonoAssignment(x))){
                  stop('"x" must contain cell-level clonotype assignment',
                       ' matrix')
              }
              summarizeClonotypes(clonoAssignment(x), by, ...)
          })







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
#'   use as this variable. Similar to the \code{group} argument for
#'   \code{\link{clonoStats}}.
#'
#' @return A list of \code{Matrix} objects providing the cell-level clonotype
#'   assignments for each unique value of \code{by} (if \code{by} denotes sample
#'   labels, each matrix in the list will contain the cells from a single
#'   sample).
#'
#' @examples
#' example(addVDJtoSCE)
#' x <- clonoStats(contigs, assignment = TRUE)
#' splitClonotypes(x, by = sce$sample)
#'
#' @import Matrix
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
#' @import Matrix
#' @export
setMethod(f = "splitClonotypes",
          signature = signature(x = "matrix"),
          definition = function(x, by){
              splitClonotypes(Matrix(x), by)
          })

#' @rdname splitClonotypes
#' @param contigs character. If \code{x} is a \code{SingleCellExperiment}, the
#'   name of the \code{SplitDataFrameList} in the \code{colData} of \code{x}
#'   containing contig information.
#' @param clonoStats character. If \code{x} is a \code{SingleCellExperiment},
#'   the name of the element in the \code{metadata} of \code{x} that contains
#'   the output of \code{clonoStats}. Must include cell-level clonotype
#'   assignments (ie. \code{assignment = TRUE}).
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
setMethod(f = "splitClonotypes",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, by, contigs = 'contigs', clonoStats = 'clonoStats'){
              if(is.null(metadata(x)[[clonoStats]])){
                  stop('No clonotype counts found.')
              }
              if(length(by) == 1){ # check SCE and contigs
                  byVar <- factor(x[[by]]) # SCE
                  if(length(byVar) == 0){
                      if(by %in% names(x[[contigs]][[1]])){
                          byVar <- factor(vapply(seq_len(ncol(x)), function(i){
                              x[[contigs]][,by][[i]][1]
                          }, FUN.VALUE = 'a')) # contigs
                      }
                  }
              }else{
                  byVar <- factor(by)
              }
              return(splitClonotypes(clonoAssignment(metadata(x)[[clonoStats]]),
                                     byVar))
          })

#' @rdname splitClonotypes
#' @export
setMethod(f = "splitClonotypes",
          signature = signature(x = "clonoStats"),
          definition = function(x, by){
              if(is.null(clonoAssignment(x))){
                  stop('"x" must contain cell-level clonotype assignment',
                       ' matrix')
              }
              splitClonotypes(clonoAssignment(x), by)
          })



#' @describeIn clonoStats-class a short summary of a \code{clonoStats}
#'   object.
#' @importFrom S4Vectors coolcat
#' @export
setMethod(f = "show",
          signature = signature(object = "clonoStats"),
          definition = function(object){
              cat('An object of class "clonoStats"\n')
              cat('clonotypes:', nrow(clonoAbundance(object)), '\n')
              cat('cells:', length(clonoGroup(object)), '\n')
              groups <- levels(clonoGroup(object))
              coolcat("groups(%d): %s\n", groups)
              cat('has assignment:', !is.null(clonoAssignment(object)))
          })

#' @describeIn clonoStats-class Get the full clonotype names
#' @export
setGeneric(name = "clonoNames",
           signature = c("object"),
           def = function(object) standardGeneric("clonoNames"))

#' @rdname clonoStats-class
#' @param object a \code{\link{clonoStats}} object or
#'   \code{SingleCellExperiment} object containing \code{clonoStats} results.
#' 
#' @export
setMethod(f = "clonoNames",
          signature = signature(object = "clonoStats"),
          definition = function(object){
              paste(object@names1, object@names2)
          })

#' @rdname clonoStats-class
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @export
setMethod(f = "clonoNames",
          signature = signature(object = "SingleCellExperiment"),
          definition = function(object){
              clonoNames(metadata(object)$clonoStats)
          })


#' @describeIn clonoStats-class Get the clonotype abundance matrix (clonotype
#'   counts per group).
#' @export
setGeneric(name = "clonoAbundance",
           signature = c("object"),
           def = function(object) standardGeneric("clonoAbundance"))

#' @rdname clonoStats-class
#' @export
setMethod(f = "clonoAbundance",
          signature = signature(object = "clonoStats"),
          definition = function(object){
              object@abundance
          })

#' @rdname clonoStats-class
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @export
setMethod(f = "clonoAbundance",
          signature = signature(object = "SingleCellExperiment"),
          definition = function(object){
              clonoAbundance(metadata(object)$clonoStats)
          })

#' @describeIn clonoStats-class Get the matrix of clonotype frequencies
#'   (singletons, doubletons, etc. per group).
#' @export
setGeneric(name = "clonoFrequency",
           signature = c("object"),
           def = function(object) standardGeneric("clonoFrequency"))

#' @rdname clonoStats-class
#' @export
setMethod(f = "clonoFrequency",
          signature = signature(object = "clonoStats"),
          definition = function(object){
              object@frequency
          })

#' @rdname clonoStats-class
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @export
setMethod(f = "clonoFrequency",
          signature = signature(object = "SingleCellExperiment"),
          definition = function(object){
              clonoFrequency(metadata(object)$clonoStats)
          })


#' @describeIn clonoStats-class Get the matrix of cell-level clonotype
#'   assignments (mapping cells to clonotypes).
#' @export
setGeneric(name = "clonoAssignment",
           signature = c("object"),
           def = function(object) standardGeneric("clonoAssignment"))

#' @rdname clonoStats-class
#' @export
setMethod(f = "clonoAssignment",
          signature = signature(object = "clonoStats"),
          definition = function(object){
              object@assignment
          })


#' @rdname clonoStats-class
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @export
setMethod(f = "clonoAssignment",
          signature = signature(object = "SingleCellExperiment"),
          definition = function(object){
              clonoAssignment(metadata(object)$clonoStats)
          })


#' @describeIn clonoStats-class Get the factor variable of group labels
#' @export
setGeneric(name = "clonoGroup",
           signature = c("object"),
           def = function(object) standardGeneric("clonoGroup"))

#' @rdname clonoStats-class
#' @export
setMethod(f = "clonoGroup",
          signature = signature(object = "clonoStats"),
          definition = function(object){
              object@group
          })


#' @rdname clonoStats-class
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @export
setMethod(f = "clonoGroup",
          signature = signature(object = "SingleCellExperiment"),
          definition = function(object){
              clonoGroup(metadata(object)$clonoStats)
          })
