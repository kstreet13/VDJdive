# quantification by Cell Ranger clonotypes
#' @title Clonotype quantification given by Cell Ranger
#' @name CRquant
#' @param ... additional arguments.
#' @export
setGeneric(name = "CRquant",
           signature = "x",
           def = function(x, ...) standardGeneric("CRquant"))

#' @rdname CRquant
#'
#' @description Use clonotype assignments produced by Cell Ranger. These are
#'   assigned separately for each sample, so clonotypes generally cannot be
#'   compared across samples.
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
#'   \code{ncol(x)}) indicating the sample of origin.
#' @param type The type of VDJ data (\code{"TCR"} or \code{"BCR"}). If
#'   \code{NULL}, this is determined by the most prevalent \code{chain} types in
#'   \code{x}.
#'   
#' @details This quantification method uses the clonotype labels provided by the
#'   Cell Ranger V(D)J pipeline. This will provide a label for every cell that
#'   has V(D)J information, but these labels are not directly comparable across
#'   samples.
#'
#' @return Creates a sparse matrix (\code{dgRMatrix}) of cell-level clonotype
#'   assignments (cells-by-clonotypes). If \code{x} is a
#'   \code{SingleCellExperiment}, this matrix is added to the \code{colData}
#'   under the name \code{clono}.
#'
#' @examples
#' data('contigs')
#' counts <- CRquant(contigs)
#'
#' @importClassesFrom IRanges SplitDataFrameList
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
setMethod(f = "CRquant",
          signature = signature(x = "SingleCellExperiment"),
          definition = function(x, TCRcol = 'contigs',
                                sample = NULL, type = NULL){
              sce <- x
              if(!is.null(sample)){
                  if(length(sample) == 1){
                      stopifnot(sample %in% names(colData(sce)))
                      sampVar <- factor(sce[[sample]])
                  }else{
                      stopifnot(length(sample) == ncol(sce))
                      sampVar <- factor(sample)
                  }
                  # calculate cells x clonotypes matrix
                  clono <- CRquant(sce[[TCRcol]], sample = sampVar, type = type)
              }else{
                  # calculate cells x clonotypes matrix
                  clono <- CRquant(sce[[TCRcol]])
              }

              # update sce
              colData(sce)$clono <- clono
              return(sce)
          })

#' @rdname CRquant
#' @importClassesFrom IRanges SplitDataFrameList
#' @importFrom Matrix Matrix colSums
#' @importClassesFrom Matrix dgRMatrix
#' @export
setMethod(f = "CRquant",
          signature = signature(x = "SplitDataFrameList"),
          definition = function(x, sample = 'sample', type = NULL){

              contigs <- x
              if(is.null(type)){
                  chn <- unlist(contigs[,'chain'])
                  if(sum(chn %in% c('IGH','IGL','IGK')) > 
                     sum(chn %in% c('TRA','TRB','TRD','TRG'))){
                      type <- 'BCR'
                  }else{
                      type <- 'TCR'
                  }
              }

              if(!is.null(sample)){
                  if(length(sample) == 1){
                      stopifnot(sample %in% colnames(contigs[[1]]))
                      sampVar <- factor(sapply(contigs[,sample], function(x){ x[1] }))
                  }else{
                      stopifnot(length(sample) == length(contigs))
                      sampVar <- factor(sample)
                  }
                  clono.list <- lapply(levels(sampVar), function(lv){
                      clono <- CRquant(contigs[which(sampVar == lv)],
                                         sample = NULL, type = type)
                      colnames(clono) <- paste0(lv, '_', colnames(clono))
                      return(clono)
                  })
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
                  all.clonotypes <- unique(unlist(lapply(clono.list, colnames)))
                  clono.list <- lapply(clono.list, function(mat){
                      missing <- all.clonotypes[! all.clonotypes %in% colnames(mat)]
                      zeros <- Matrix(0, ncol = length(missing), nrow = nrow(mat),
                                      sparse = TRUE)
                      colnames(zeros) <- missing
                      return(cbind(mat, zeros))
                  })
                  clono <- do.call(rbind, clono.list)
                  clono <- clono[names(contigs), ]
                  return(clono)
              }

              # remove unproductive and 'Multi' contigs (for now?)
              contigs <- contigs[contigs[,'productive']]
              if(type == 'TCR'){
                  type1 <- 'TRA'
                  type2 <- 'TRB'
              }
              if(type == 'BCR'){
                  type1 <- 'IGH'
                  type2 <- c('IGL','IGK')
                  # prepend "IGK" or "IGL" to distinguish?
                  cellbarcodes <- names(contigs)
                  contigs <- unlist(contigs)
                  t2ind <- which(contigs$chain %in% type2)
                  contigs$cdr3[t2ind] <- paste(contigs$chain[t2ind], 
                                               contigs$cdr3[t2ind], sep = '-')
                  contigs <- split(DataFrame(contigs), 
                                   factor(contigs$barcode, cellbarcodes))
              }
              contigs <- contigs[contigs[,'chain'] %in% c(type1, type2)]
              

              all_clonotypes <- unique(unlist(contigs[,'raw_clonotype_id']))

              poss.indices <- sapply(contigs, function(x){
                  which(all_clonotypes == x$raw_clonotype_id[1])
              })


              # build full cells-by-clonotypes matrix
              clonoJ <- as.integer(unlist(poss.indices)-1)
              clonoP <- as.integer(c(0,cumsum(lengths(poss.indices))))
              clonoX <- rep(1, sum(lengths(poss.indices)))
              clono <- new('dgRMatrix', j = clonoJ, p = clonoP, x = clonoX,
                           Dim = as.integer(c(length(contigs), length(all_clonotypes))))
              colnames(clono) <- all_clonotypes
              clono <- clono[, which(colSums(clono) > 0), drop = FALSE]
              rownames(clono) <- names(contigs)

              return(clono)
          })

