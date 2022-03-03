# quantification for cells that can be uniquely assigned
# (ie. have exactly one alpha and one beta chain)

#' @title Clonotype quantification for uniquely assignable cells
#' @name uniqueQuant
#' @param ... additional arguments.
#' @export
setGeneric(name = "uniqueQuant",
           signature = "x",
           def = function(x, ...) standardGeneric("uniqueQuant"))

#' @rdname uniqueQuant
#'
#' @description Assign cells that can be uniquely mapped to a single clonotype.
#'   This results in fewer cells being used compared to other quantification
#'   methods, but all assignments are unambiguous.
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
#' @details This quantification method defines a clonotype as a pair of specific
#'   chains (alpha and beta for T cells, heavy and light for B cells) and only
#'   assigns cells that are a perfect match for their clonotype, leaving others
#'   unassigned, even if they contain V(D)J information.
#'
#' @return Creates a sparse matrix (\code{dgRMatrix}) of cell-level clonotype
#'   assignments (cells-by-clonotypes). If \code{x} is a
#'   \code{SingleCellExperiment}, this matrix is added to the \code{colData}
#'   under the name \code{clono}.
#'
#' @examples
#' data('contigs')
#' counts <- uniqueQuant(contigs)
#'
#' @importClassesFrom IRanges SplitDataFrameList
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
setMethod(f = "uniqueQuant",
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
                  clono <- uniqueQuant(sce[[TCRcol]], 
                                       sample = sampVar, type = type)
              }else{
                  # calculate cells x clonotypes matrix
                  clono <- uniqueQuant(sce[[TCRcol]], type = type)
              }

              # update sce
              colData(sce)$clono <- clono
              return(sce)
          })

#' @rdname uniqueQuant
#' @importClassesFrom IRanges SplitDataFrameList
#' @importFrom Matrix Matrix colSums
#' @importClassesFrom Matrix dgRMatrix
#' @export
setMethod(f = "uniqueQuant",
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
                      uniqueQuant(contigs[which(sampVar == lv)],
                                  sample = NULL, type = type)
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

              
              # explore possible numbers of alpha and beta chains
              nAlpha <- sum(contigs[,"chain"] %in% type1)
              nBeta <- sum(contigs[,"chain"] %in% type2)
              
              # find all unique alpha chains
              all.alphas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']%in%type1]))
              if(length(all.alphas)==0){
                  all.alphas <- 'unknown'
              }
              # find all unique beta chains
              all.betas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']%in%type2]))
              if(length(all.betas)==0){
                  all.betas <- 'unknown'
              }
              
              # initialize counts matrix (#alpha-by-#beta)
              counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), sparse = TRUE)
              rownames(counts) <- all.alphas
              colnames(counts) <- all.betas

              # useful indices
              ind.unique <- which(nAlpha == 1 & nBeta == 1)

              # use which alpha sequence and which beta sequence each contig represents to calculate its (possible) index in the clonotype matrix
              wa <- match(contigs[,'cdr3'][contigs[,'chain']%in%type1], all.alphas)
              wb <- match(contigs[,'cdr3'][contigs[,'chain']%in%type2], all.betas)
              # this takes care of the 1-alpha + 1-beta indices without looping
              poss.indices <- as.list(length(all.alphas)*(wb-1L) + wa)
              # others are incorrect, though
              poss.indices[-ind.unique] <- list(integer(0))

              # step 1 (the only step): assign cells with 1 alpha, 1 beta
              ############################################
              temp <- contigs[ind.unique]
              uniquecounts <- unclass(table(unlist(temp[temp[,'chain']%in%type1,'cdr3']),
                                            unlist(temp[temp[,'chain']%in%type2,'cdr3'])))
              if(length(temp) > 0){
                  counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
              }
              counts <- as.numeric(counts)


              # build full cells-by-clonotypes matrix
              clonoJ <- as.integer(unlist(poss.indices)-1)
              clonoP <- as.integer(c(0,cumsum(lengths(poss.indices))))
              clonoX <- rep(1, sum(lengths(poss.indices)))
              clono <- new('dgRMatrix', j = clonoJ, p = clonoP, x = clonoX,
                           Dim = as.integer(c(length(contigs), length(counts))))
              colnames(clono) <- as.character(outer(all.alphas, all.betas,
                                                    FUN = paste))
              clono <- clono[, which(colSums(clono) > 0), drop = FALSE]
              rownames(clono) <- names(contigs)

              return(clono)
          })



