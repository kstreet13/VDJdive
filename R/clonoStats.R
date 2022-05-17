#' @include clonoStats_class.R
NULL

#' @title Assign cell-level clonotypes and calculate abundances
#' @param ... additional arguments.
#' @name clonoStats
#' @export
setGeneric(name = "clonoStats",
           signature = "x",
           def = function(x, ...) standardGeneric("clonoStats"))

#' @rdname clonoStats
#'
#' @description Assign clonotypes to cells and produce two summary tables: the
#'   \code{clonotypes x samples} table of abundances and the \code{counts x
#'   samples} table of clonotype frequencies.
#'
#' @param x A \code{SplitDataFrameList} object containing V(D)J contig
#'   information, split by cell barcodes, as created by \code{readVDJcontigs}.
#'   Alternatively, a \code{SingleCellExperiment} object with such a
#'   \code{SplitDataFrameList} in the \code{colData}, as created by
#'   \code{addVDJtoSCE}.
#' @param TCRcol The name of the column in the \code{colData} of \code{x} that
#'   contains the VDJ contig data (only applies if \code{x} is a
#'   \code{SingleCellExperiment}).
#' @param sample The name of the column in \code{x} (or the \code{colData} of
#'   \code{x}, for \code{SingleCellExperiment} objects) that stores each cell's
#'   sample of origin. Alternatively, a vector of length equal to \code{x} (or
#'   \code{ncol(x)}) indicating the sample of origin. Providing this information
#'   can dramatically speed up computation and avoid unwanted cross-talk between
#'   samples.
#' @param type The type of VDJ data (\code{"TCR"} or \code{"BCR"}). If
#'   \code{NULL}, this is determined by the most prevalent \code{chain} types in
#'   \code{x}.
#' @param assignment Logical, whether or not to return the full \code{nCells x
#'   nClonotypes} sparse matrix of clonotype assignments (default =
#'   \code{FALSE})
#' @param method Which method to use for assigning cell-level clonotypes.
#'   Options are \code{"EM"} (default), \code{"unique"}, \code{"CellRanger"}, or
#'   any \code{chain} type in \code{x}. See Details.
#' @param lang Indicates which implementation of certain methods to use. The EM
#'   algorithm is implemented in both pure R (\code{'r'}) and mixed R and Python
#'   (\code{'python'}, default) versions. Similarly, clonotype summarization is
#'   implemented in two ways, which can impact speed, regardless of choice of
#'   \code{method}.
#' @param thresh Numeric threshold for convergence of the EM algorithm.
#'   Indicates the maximum allowable deviation in a count between updates. Only
#'   used if \code{method = "EM"}.
#' @param iter.max Maximum number of iterations for the EM algorithm. Only used
#'   if \code{method = "EM"}.
#'
#' @details Assign each cell (with at least one V(D)J contig) to its most
#'   likely clonotype with the EM algorithm. For ambiguous cells, this leads to
#'   proportional (non-integer) assignment across multiple possible clonotypes.

#' @details This quantification method defines a clonotype as a pair of specific
#'   chains (alpha and beta for T cells, heavy and light for B cells) and
#'   attempts to assign each cell to its most likely clonotype. Unlike other
#'   quantification methods, this can lead to non-integer counts for cells with
#'   ambiguous information (ie. only an alpha chain, or two alphas and one beta
#'   chain).
#'
#' @details We highly recommend providing information on each cell's sample of
#'   origin, as this can speed up computation and provide more accurate results.
#'   Because the EM algorithm shares information across cells, splitting by
#'   sample can improve accuracy by removing extraneous clonotypes from the set
#'   of possibilities for a particular cell.
#'
#' @return Creates a sparse matrix (\code{dgRMatrix}) of cell-level clonotype
#'   assignments (cells-by-clonotypes). If \code{x} is a
#'   \code{SingleCellExperiment}, this matrix is added to the \code{colData}
#'   under the name \code{clono}.
#'
#' @examples
#' data('contigs')
#' clonoStats(contigs)
#'
#' @import IRanges
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @importFrom Matrix Matrix colSums
#' @importClassesFrom Matrix dgRMatrix
#'
#' @export
setMethod(f = "clonoStats",
          signature = signature(x = "SplitDataFrameList"),
          definition = function(x, sample = 'sample', type = NULL,
                                assignment = FALSE, 
                                method = 'EM', 
                                lang = c('python','r'),
                                thresh = .01, iter.max = 1000){
              contigs <- x
              method <- match.arg(method, 
                                  choices = c('EM','unique','CellRanger',
                                              unique(unlist(x[,'chain']))))
              lang <- match.arg(lang)
              if(is.null(type)){
                  chn <- unlist(contigs[,'chain'])
                  if(sum(chn %in% c('IGH','IGL','IGK')) > 
                     sum(chn %in% c('TRA','TRB','TRD','TRG'))){
                      type <- 'BCR'
                  }else{
                      type <- 'TCR'
                  }
              }
              
              if(is.null(sample)){
                  sampVar <- factor(rep('sample1', length(contigs)))
              }else{
                  if(length(sample) == 1){
                      stopifnot(sample %in% colnames(contigs[[1]]))
                      sampVar <- factor(vapply(seq_along(contigs), function(i){
                          contigs[[i]][,sample][1]
                      }, FUN.VALUE = 'a'))
                  }else{
                      if(!is.factor(sample)){
                          sampVar <- factor(sample)
                      }else{ # keep samples with no TCR data
                          # (factor() drops empty levels)
                          sampVar <- sample
                      }
                      stopifnot(length(sampVar) == length(contigs))
                  }
              }
              
              # select appropriate method and quantify each sample individually
              #################################################################
              if(method == 'EM'){
                  clono.list <- lapply(levels(sampVar), function(lv){
                      .EM_sample(contigs[which(sampVar == lv)],
                                 type = type, lang = lang,
                                 thresh = thresh, iter.max = iter.max)
                  })
              }else if(method == 'CellRanger'){
                  clono.list <- lapply(levels(sampVar), function(lv){
                      .CR_sample(contigs[which(sampVar == lv)],
                                 type = type)
                  })
              }else if(method == 'unique'){
                  clono.list <- lapply(levels(sampVar), function(lv){
                      .UNIQ_sample(contigs[which(sampVar == lv)],
                                   type = type)
                  })
              }else{
                  clono.list <- lapply(levels(sampVar), function(lv){
                      .CHN_sample(contigs[which(sampVar == lv)], 
                                  method = method)
                  })
              }
              if(any(is.na(sampVar))){
                  # indicates elements of length 0
                  clono.list[[length(clono.list)+1]] <-
                      Matrix(0, nrow = sum(is.na(sampVar)),
                             ncol = ncol(clono.list[[1]]), sparse = TRUE)
                  rownames(clono.list[[length(clono.list)]]) <-
                      names(contigs)[is.na(sampVar)]
                  colnames(clono.list[[length(clono.list)]]) <-
                      colnames(clono.list[[1]])
              }
              # supplement each matrix with 0s so they have same number of cols
              all.clonotypes <- unique(unlist(lapply(clono.list, colnames)))
              if(length(all.clonotypes) > 0){ # can't index by empty set
                  clono.list <- lapply(clono.list, function(mat){
                      missing <- all.clonotypes[! all.clonotypes %in% 
                                                    colnames(mat)]
                      zeros <- Matrix(0, ncol = length(missing), 
                                      nrow = nrow(mat), sparse = TRUE)
                      colnames(zeros) <- missing
                      return(cbind(mat, zeros)[,all.clonotypes])
                  })
              }
              # combine matrices
              clono <- do.call(rbind, clono.list)
              clono <- clono[names(contigs), ]
              
              # make names factors
              if(is.null(colnames(clono))){
                  nf1 <- nf2 <- factor()
              }else{
                  s <- strsplit(colnames(clono), split=' ')
                  ind <- which(lengths(s) == 1)
                  for(i in ind){
                      s[[i]] <- c(s[[i]], '')
                  }
                  alphas <- unlist(s)[seq(1,2*ncol(clono), by=2)]
                  betas <- unlist(s)[seq(2,2*ncol(clono), by=2)]
                  nf1 <- factor(alphas)
                  nf2 <- factor(betas)
                  colnames(clono) <- NULL
              }
              
              # summarize clonotype assignments
              t1 <- summarizeClonotypes(clono, sampVar, mode = 'sum')
              rownames(t1) <- NULL
              t2 <- summarizeClonotypes(clono, sampVar, mode = 'tab', 
                                        lang = lang)
              if(assignment){
                  return(new('clonoStats', abundance = t1, frequency = t2,
                             names1 = nf1, names2 = nf2, group = sampVar,
                             assignment = clono))
              }else{
                  return(new('clonoStats', abundance = t1, frequency = t2,
                             names1 = nf1, names2 = nf2, group = sampVar))
              }
          })

#' @rdname clonoStats
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata<-
#' @export
setMethod(f = "clonoStats",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, TCRcol = 'contigs', sample = NULL, ...){
              sce <- x
              if(!is.null(sample)){
                  if(length(sample) == 1){
                      stopifnot(sample %in% names(colData(sce)))
                      sampVar <- factor(sce[[sample]])
                  }else{
                      stopifnot(length(sample) == ncol(sce))
                      sampVar <- factor(sample)
                  }
                  cs <- clonoStats(sce[[TCRcol]], 
                                   sample = sampVar, ...)
              }else{
                  cs <- clonoStats(sce[[TCRcol]], ...)
              }
              
              # update sce
              # if(!is.null(clono@assignment)){
              #     colData(sce)$assignment <- clono@assignment
              #     clono@assignment <- NULL
              # }
              metadata(sce)$clonoStats <- cs
              
              return(sce)
          })

#' @rdname clonoStats
#' @export
setMethod(f = "clonoStats",
          signature = signature(x = "clonoStats"),
          definition = function(x, sample = NULL, lang = c('python','r')){
              # remake clonoStats object with new group variable
              lang <- match.arg(lang)
              if(is.null(x@assignment)){
                  stop('"x" must contain cell-level clonotype assignments')
              }
              if(is.null(sample)){
                  sampVar <- factor(rep(1, length(x@group)))
              }else{
                  stopifnot(length(sample) == length(x@group))
                  sampVar <- factor(sample)
              }
              stopifnot(length(sampVar) == nrow(x@assignment))
              # summaries
              t1 <- summarizeClonotypes(x@assignment, sampVar, mode = 'sum')
              rownames(t1) <- NULL
              t2 <- summarizeClonotypes(x@assignment, sampVar, mode = 'tab', 
                                        lang = lang)
              # keep assignment matrix
              return(new('clonoStats', abundance = t1, frequency = t2,
                         names1 = x@names1, names2 = x@names2, group = sampVar,
                         assignment = x@assignment))
          })





