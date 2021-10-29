

#' @param sce SingleCellExperiment object with TCR data
#' @param TCRcol element of the colData of sce that contains the TCR data
#' @param sample optional variable for splitting quantification by sample. This can speed up computation and avoid crosstalk between samples.
#' @param thresh threshold for convergence in assigning clonotype counts
#' @param iter.max maximum number of iterations for the EM algorithm.
require(Matrix)
require(SingleCellExperiment)
require(IRanges)
reticulate::source_python('src/updateCounts.py')

# next step: parallelize it

setGeneric(name = "EMquant",
           signature = "x",
           def = function(x, ...) standardGeneric("EMquant"))

setMethod(f = "EMquant",
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
                  # calculate cells x clonotypes matrix
                  clono <- EMquant(sce[[TCRcol]], sample = sampVar, ...)
              }else{
                  # calculate cells x clonotypes matrix
                  clono <- EMquant(sce[[TCRcol]], ...)
              }
              
              # update sce
              colData(sce)$clono <- clono
              return(sce)
          })

setMethod(f = "EMquant",
          signature = signature(x = "SplitDataFrameList"),
          definition = function(x, sample = 'sample', 
                                method = c('python','r'),
                                thresh = .01, iter.max = 1000){
              contigs <- x
              
              method <- match.arg(method)
              if(!is.null(sample)){
                  if(length(sample) == 1){
                      stopifnot(sample %in% colnames(contigs[[1]]))
                      sampVar <- factor(sapply(contigs[,sample], function(x){ x[1] }))
                  }else{
                      stopifnot(length(sample) == length(contigs))
                      sampVar <- factor(sample)
                  }
                  clono.list <- lapply(levels(sampVar), function(lv){
                      EMquant(contigs[which(sampVar == lv)],
                              sample = NULL, method = method,
                              thresh = thresh, iter.max = iter.max)
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
              contigs <- contigs[contigs[,'productive']=='True']
              contigs <- contigs[contigs[,'chain'] %in% c('TRA','TRB')]
              
              # explore possible numbers of alpha and beta chains
              nAlpha <- sum(contigs[,"chain"] == 'TRA')
              nBeta <- sum(contigs[,"chain"] == 'TRB')
              
              # find all unique alpha chains
              all.alphas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRA']))
              if(length(all.alphas)==0){
                  all.alphas <- 'unknown'
              }
              # find all unique beta chains
              all.betas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRB']))
              if(length(all.betas)==0){
                  all.betas <- 'unknown'
              }
              
              # initialize counts matrix (#alpha-by-#beta)
              require(Matrix)
              counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), sparse = TRUE)
              rownames(counts) <- all.alphas
              colnames(counts) <- all.betas
              
              # useful indices
              ind.unique <- which(nAlpha == 1 & nBeta == 1)
              ind.multiple <- which(nAlpha > 0 & nBeta > 0 & nAlpha+nBeta > 2)
              ind.noAlpha <- which(nAlpha == 0 & nBeta > 0)
              ind.noBeta <- which(nBeta == 0 & nAlpha > 0)
              ind.ambiguous <- sort(c(ind.multiple, ind.noAlpha, ind.noBeta))
              
              # use which alpha sequence and which beta sequence each contig represents to calculate its (possible) index in the clonotype matrix
              wa <- match(contigs[,'cdr3'][contigs[,'chain']=='TRA'], all.alphas)
              wb <- match(contigs[,'cdr3'][contigs[,'chain']=='TRB'], all.betas)
              # this takes care of the 1-alpha + 1-beta indices without looping
              poss.indices <- as.list(length(all.alphas)*(wb-1L) + wa)
              # others are incorrect, though
              poss.indices[ind.multiple] <- lapply(ind.multiple, function(i){
                  return(sort(as.integer(sapply(wa[[i]], function(a){
                      sapply(wb[[i]], function(b){
                          length(all.alphas)*(b-1) + a
                      })
                  }))))
              })
              poss.indices[ind.noAlpha] <- lapply(ind.noAlpha, function(i){
                  a.i <- seq_along(all.alphas)
                  return(sort(as.integer(sapply(wb[[i]], function(b){
                      length(all.alphas)*(b-1) + a.i
                  }))))
              })
              poss.indices[ind.noBeta] <- lapply(ind.noBeta, function(i){
                  b.i <- seq_along(all.betas)
                  return(sort(as.integer(sapply(wa[[i]], function(a){
                      sort(length(all.alphas)*(b.i-1) + a)
                  }))))
              })
              
              # step 1: assign cells with 1 alpha, 1 beta
              ############################################
              temp <- contigs[ind.unique]
              uniquecounts <- unclass(table(unlist(temp[temp[,'chain']=='TRA','cdr3']), 
                                            unlist(temp[temp[,'chain']=='TRB','cdr3'])))
              if(length(temp) > 0){
                  counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
              }
              uniquecounts <- counts <- as.numeric(counts)
              
              # step 2: assign (multi-mapped) cells (equally across all possibilities)
              ############################################################
              temp <- contigs[ind.ambiguous]
              if(length(temp) > 0){
                  t.indices <- poss.indices[ind.ambiguous]
                  for(i in seq_along(temp)){
                      counts[t.indices[[i]]] <- 
                          counts[t.indices[[i]]] + 1/length(t.indices[[i]])
                  }
                  counts.old <- counts
              }
              
              # repeat 2 (proportional to previous counts)
              #################
              if(method == 'python'){
                  t.indices <- poss.indices[ind.ambiguous]
                  names(t.indices) <- NULL # so python takes it as a list of lists, not a dict
                  # force python to take length-1 vectors as lists
                  l1.idx <- which(lengths(t.indices) == 1)
                  for(ii in l1.idx){
                      t.indices[[ii]] <- list(as.integer(t.indices[[ii]]))
                  }
                  
                  # iteration handled by python
                  counts <- TCR_EM_counts(uniquecounts, counts.old, t.indices, thresh, iter.max)
              }else{
                  working <- TRUE
                  iters <- 0
                  temp <- contigs[ind.ambiguous]
                  t.indices <- poss.indices[ind.ambiguous]
                  while(working){
                      iters <- iters + 1
                      counts <- uniquecounts
                      
                      # step 2 (proportional to counts.old)
                      #########
                      for(i in seq_along(temp)){
                          counts[t.indices[[i]]] <- 
                              counts[t.indices[[i]]] + 
                              counts.old[t.indices[[i]]] / 
                              sum(counts.old[t.indices[[i]]])
                      }
                      
                      # check for convergence
                      diff <- max(abs(counts-counts.old))
                      #print(diff[length(diff)])
                      if(diff < thresh | iters >= iter.max){
                          working <- FALSE
                      }else{
                          counts.old <- counts
                      }
                  }
              }
              
              # build full cells-by-clonotypes matrix
              clonoJ <- as.integer(unlist(poss.indices)-1)
              clonoP = as.integer(c(0,cumsum(lengths(poss.indices))))
              clonoX <- unlist(sapply(which(lengths(poss.indices) > 0), function(i){
                  counts[poss.indices[[i]]] / 
                      sum(counts[poss.indices[[i]]])
              }))
              clono <- new('dgRMatrix', j = clonoJ, p = clonoP, x = clonoX,
                           Dim = as.integer(c(length(contigs), length(counts))))
              colnames(clono) <- as.character(outer(all.alphas, all.betas, 
                                                    FUN = paste))
              clono <- clono[, colSums(clono) > 0]
              rownames(clono) <- names(contigs)
              
              return(clono)
          })


