#' @include clonoStats_class.R clonoStats_helpers.R
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
#' @description Assign clonotype labels to cells and produce two summary tables:
#'   the \code{clonotypes x samples} table of abundances and the \code{counts x
#'   samples} table of clonotype frequencies.
#'
#' @param x A \code{SplitDataFrameList} object containing V(D)J contig
#'   information, split by cell barcodes, as created by \code{readVDJcontigs}.
#'   Alternatively, a \code{SingleCellExperiment} object with such a
#'   \code{SplitDataFrameList} in the \code{colData}, as created by
#'   \code{addVDJtoSCE}.
#' @param contigs character. When \code{x} is a \code{SingleCellExperiment},
#'   this is the name of the column in the \code{colData} of \code{x} that
#'   contains the VDJ contig data.
#' @param group character. The name of the column in \code{x} (or the
#'   \code{colData} of \code{x}, for \code{SingleCellExperiment} objects) that
#'   stores each cell's group identity, typically either its sample of origin or
#'   cluster label. Alternatively, a vector of length equal to \code{x} (or
#'   \code{ncol(x)}) indicating the group identity. Providing this information
#'   can dramatically speed up computation. When running \code{clonoStats} for
#'   the first time on a dataset, we highly recommend setting the group identity
#'   to sample of origin to avoid unwanted cross-talk between samples.
#' @param type character. The type of VDJ data (one of \code{"TCR"} or
#'   \code{"BCR"}). If \code{NULL}, this is determined by the most prevalent
#'   \code{chain} types in \code{x}.
#' @param assignment logical. Whether or not to return the full \code{nCells x
#'   nClonotypes} sparse matrix of clonotype assignments (default =
#'   \code{FALSE})
#' @param method character. Which method to use for assigning cell-level
#'   clonotypes. Options are \code{"EM"} (default), \code{"unique"}, or
#'   \code{"CellRanger"}. Alternatively, this may be the name of a numeric
#'   column of the contig data or any \code{chain} type contained therein. See
#'   Details.
#' @param lang character. Indicates which implementation of certain methods to
#'   use. The EM algorithm is implemented in both pure R (\code{'r'}) and mixed
#'   R and Python (\code{'python'}, default) versions. Similarly, clonotype
#'   summarization is implemented in two ways, which can impact speed,
#'   regardless of choice of \code{method}.
#' @param thresh Numeric threshold for convergence of the EM algorithm.
#'   Indicates the maximum allowable deviation in a count between updates. Only
#'   used if \code{method = "EM"}.
#' @param iter.max Maximum number of iterations for the EM algorithm. Only used
#'   if \code{method = "EM"}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the
#'   parallel backend for distributed clonotype assignment operations (split by
#'   \code{group}). Default is \code{BiocParallel::SerialParam()}.
#'
#' @details Assign cells (with at least one V(D)J contig) to clonotypes and
#'   produce summary tables that can be used for downstream analysis. Clonotype
#'   assignment can be handled in multiple ways depending on the choice of
#'   \code{"method"}:
#'   \itemize{
#'   \item{\code{"EM"}: }{Cells are assigned probabilistically to their most
#'   likely clonotype(s) with the Expectation-Maximization (EM) algorithm. For
#'   ambiguous cells, this leads to proportional (non-integer) assignment across
#'   multiple clonotypes and a frequency table of (non-integer) expected
#'   counts.}
#'   \item{\code{"unique"}: }{Cells are assigned a clonotype if (and only if)
#'   they can be uniquely assigned a single clonotype. For a T cell, this means
#'   having exactly one alpha chain and one beta chain.}
#'   \item{\code{"CellRanger"}: }{Clonotype labels are taken from contig data
#'   and matched across samples.}
#'   \item{\code{column name in contig data}: }{Similar to \code{"unique"}, but
#'   additionally, cells with multiples of a particular chain are assigned a
#'   "dominant" clonotype based on which contig has the higher value in this
#'   column (typical choices being \code{"umis"} or \code{"reads"}).}
#'   \item{\code{type of chain in contig data}: }{Clonotypes are based entirely
#'   on this type of chain (eg. \code{"TRA"} or \code{"TRB"}) and cells may be
#'   assigned to multiple clonotypes, if multiples of that chain are present.} }
#'   
#' @details The \code{"EM"}, \code{"unique"}, and UMI/read-based quantification
#'   methods all define a clonotype as a pair of specific chains (alpha and beta
#'   for T cells, heavy and light for B cells). Unlike other methods, the EM
#'   algorithm assigns clonotypes probabilistically, which can lead to
#'   non-integer counts for cells with ambiguous information (ie. only an alpha
#'   chain, or two alphas and one beta chain).
#'
#' @details We highly recommend providing information on each cell's sample of
#'   origin, as this can speed up computation and provide more accurate results.
#'   This is particularly important for the EM algorithm, which shares
#'   information across cells in the same group, so splitting by sample can
#'   improve accuracy by removing extraneous clonotypes from the set of
#'   possibilities for a particular cell.
#'
#' @return Returns an object of class \code{clonoStats}, containing group-level
#'   clonotype summaries. May optionally include a sparse matrix of cell-level
#'   assignment information, if \code{assignment = TRUE}. If \code{x} is a
#'   \code{SingleCellExperiment} object, this output is added to the metadata.
#'
#' @seealso \code{\linkS4class{clonoStats}}
#'
#' @examples
#' data('contigs')
#' clonoStats(contigs)
#'
#' @importClassesFrom IRanges SplitDataFrameList
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @importFrom Matrix Matrix colSums
#' @importClassesFrom Matrix dgRMatrix
#'
#' @export
setMethod(f = "clonoStats",
          signature = signature(x = "SplitDataFrameList"),
          definition = function(x, group = 'sample', type = NULL,
                                assignment = FALSE, 
                                method = 'EM', 
                                lang = c('python','r'),
                                thresh = .01, iter.max = 1000,
                                BPPARAM = SerialParam()){
              contigs <- x
              method <- match.arg(method, 
                                  choices = c('EM','unique','CellRanger',
                                              unique(unlist(x[,'chain'])),
                                              commonColnames(x)))
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
              
              if(is.null(group)){
                  grpVar <- factor(rep('sample1', length(contigs)))
              }else{
                  if(length(group) == 1){
                      stopifnot(group %in% colnames(contigs[[1]]))
                      grpVar <- factor(vapply(seq_along(contigs), function(i){
                          contigs[[i]][,group][1]
                      }, FUN.VALUE = 'a'))
                  }else{
                      if(!is.factor(group)){
                          grpVar <- factor(group)
                      }else{ # keep samples with no TCR data
                          # (factor() drops empty levels)
                          grpVar <- group
                      }
                      stopifnot(length(grpVar) == length(contigs))
                  }
              }
              
              # select appropriate method and quantify each sample individually
              #################################################################
              if(method == 'EM'){
                  clono.list <- bplapply(levels(grpVar), function(lv){
                      .EM_sample(contigs[which(grpVar == lv)],
                                 type = type, lang = lang,
                                 thresh = thresh, iter.max = iter.max)
                  }, BPPARAM = BPPARAM)
              }else if(method == 'CellRanger'){
                  clono.list <- bplapply(levels(grpVar), function(lv){
                      .CR_sample(contigs[which(grpVar == lv)],
                                 type = type)
                  }, BPPARAM = BPPARAM)
              }else if(method == 'unique'){
                  clono.list <- bplapply(levels(grpVar), function(lv){
                      .UNIQ_sample(contigs[which(grpVar == lv)],
                                   type = type)
                  }, BPPARAM = BPPARAM)
              }else if(method %in% unlist(x[,'chain'])){
                  clono.list <- bplapply(levels(grpVar), function(lv){
                      .CHN_sample(contigs[which(grpVar == lv)], 
                                  method = method)
                  }, BPPARAM = BPPARAM)
              }else if(method %in% commonColnames(x)){
                  clono.list <- bplapply(levels(grpVar), function(lv){
                      .MC_sample(contigs[which(grpVar == lv)], 
                                  type = type, method = method)
                  }, BPPARAM = BPPARAM)
              }
              if(any(is.na(grpVar))){
                  # indicates elements of length 0
                  clono.list[[length(clono.list)+1]] <-
                      Matrix(0, nrow = sum(is.na(grpVar)),
                             ncol = ncol(clono.list[[1]]), sparse = TRUE)
                  rownames(clono.list[[length(clono.list)]]) <-
                      names(contigs)[is.na(grpVar)]
                  colnames(clono.list[[length(clono.list)]]) <-
                      colnames(clono.list[[1]])
              }
              # supplement each matrix with 0s so they have same number of cols
              all.clonotypes <- unique(unlist(lapply(clono.list, colnames)))
              if(length(all.clonotypes) > 0){ # can't index by empty set
                  clono.list <- bplapply(clono.list, function(mat){
                      missing <- all.clonotypes[! all.clonotypes %in% 
                                                    colnames(mat)]
                      zeros <- Matrix(0, ncol = length(missing), 
                                      nrow = nrow(mat), sparse = TRUE)
                      colnames(zeros) <- missing
                      return(cbind(mat, zeros)[,all.clonotypes])
                  }, BPPARAM = BPPARAM)
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
              t1 <- summarizeClonotypes(clono, grpVar, mode = 'sum',
                                        BPPARAM = BPPARAM)
              rownames(t1) <- NULL
              t2 <- summarizeClonotypes(clono, grpVar, mode = 'tab', 
                                        lang = lang, BPPARAM = BPPARAM)
              if(assignment){
                  return(new('clonoStats', abundance = t1, frequency = t2,
                             names1 = nf1, names2 = nf2, group = grpVar,
                             assignment = clono))
              }else{
                  return(new('clonoStats', abundance = t1, frequency = t2,
                             names1 = nf1, names2 = nf2, group = grpVar))
              }
          })

#' @rdname clonoStats
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata<-
#' @export
setMethod(f = "clonoStats",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, contigs = 'contigs', group = 'sample', ...){
              sce <- x
              if(!is.null(group)){
                  if(length(group) == 1){ # check SCE and contigs
                      grpVar <- factor(sce[[group]]) # SCE
                      if(length(grpVar) == 0){
                          if(group %in% names(x[[contigs]][[1]])){
                              grpVar <- factor(vapply(seq_len(ncol(x)),
                                                      function(i){
                                  x[[contigs]][,group][[i]][1]
                              }, FUN.VALUE = 'a')) # contigs
                          }
                      }
                  }else{
                      stopifnot(length(group) == ncol(sce))
                      grpVar <- factor(group)
                  }
                  cs <- clonoStats(sce[[contigs]], 
                                   group = grpVar, ...)
              }else{
                  cs <- clonoStats(sce[[contigs]], ...)
              }
              
              metadata(sce)$clonoStats <- cs
              
              return(sce)
          })

#' @rdname clonoStats
#' @export
setMethod(f = "clonoStats",
          signature = signature(x = "clonoStats"),
          definition = function(x, group = NULL, lang = c('python','r')){
              # remake clonoStats object with new group variable
              lang <- match.arg(lang)
              if(is.null(clonoAssignment(x))){
                  stop('"x" must contain cell-level clonotype assignments')
              }
              if(is.null(group)){
                  # grpVar <- clonoGroup(x)
                  return(x)
              }else{
                  stopifnot(length(group) == length(clonoGroup(x)))
                  grpVar <- factor(group)
              }
              stopifnot(length(grpVar) == nrow(clonoAssignment(x)))
              # summaries
              t1 <- summarizeClonotypes(clonoAssignment(x), grpVar, 
                                        mode = 'sum')
              rownames(t1) <- NULL
              t2 <- summarizeClonotypes(clonoAssignment(x), grpVar, 
                                        mode = 'tab', lang = lang)
              # keep assignment matrix
              return(new('clonoStats', abundance = t1, frequency = t2,
                         names1 = x@names1, names2 = x@names2, group = grpVar,
                         assignment = clonoAssignment(x)))
          })





