# other quantification methods


# quantification for cells that can be uniquely assigned
# (ie. have exactly one alpha and one beta chain)
setGeneric(name = "uniqueQuant",
           signature = "x",
           def = function(x, ...) standardGeneric("uniqueQuant"))

setMethod(f = "uniqueQuant",
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
                  clono <- uniqueQuant(sce[[TCRcol]], sample = sampVar, ...)
              }else{
                  # calculate cells x clonotypes matrix
                  clono <- uniqueQuant(sce[[TCRcol]], ...)
              }
              
              # update sce
              colData(sce)$clono <- clono
              return(sce)
          })

setMethod(f = "uniqueQuant",
          signature = signature(x = "SplitDataFrameList"),
          definition = function(x, sample = 'sample'){
              
              contigs <- x
              
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
                                  sample = NULL)
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
              
              # use which alpha sequence and which beta sequence each contig represents to calculate its (possible) index in the clonotype matrix
              wa <- match(contigs[,'cdr3'][contigs[,'chain']=='TRA'], all.alphas)
              wb <- match(contigs[,'cdr3'][contigs[,'chain']=='TRB'], all.betas)
              # this takes care of the 1-alpha + 1-beta indices without looping
              poss.indices <- as.list(length(all.alphas)*(wb-1L) + wa)
              # others are incorrect, though
              poss.indices[-ind.unique] <- list(integer(0))

              # step 1 (the only step): assign cells with 1 alpha, 1 beta
              ############################################
              temp <- contigs[ind.unique]
              uniquecounts <- unclass(table(unlist(temp[temp[,'chain']=='TRA','cdr3']), 
                                            unlist(temp[temp[,'chain']=='TRB','cdr3'])))
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
              clono <- clono[, colSums(clono) > 0]
              rownames(clono) <- names(contigs)
              
              return(clono)
          })




# quantification by whatever Cell Ranger is doing
setGeneric(name = "tenXquant",
           signature = "x",
           def = function(x, ...) standardGeneric("tenXquant"))

setMethod(f = "tenXquant",
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
                  clono <- tenXquant(sce[[TCRcol]], sample = sampVar, ...)
              }else{
                  # calculate cells x clonotypes matrix
                  clono <- tenXquant(sce[[TCRcol]], ...)
              }
              
              # update sce
              colData(sce)$clono <- clono
              return(sce)
          })

setMethod(f = "tenXquant",
          signature = signature(x = "SplitDataFrameList"),
          definition = function(x, sample = 'sample'){
              
              contigs <- x
              
              if(!is.null(sample)){
                  if(length(sample) == 1){
                      stopifnot(sample %in% colnames(contigs[[1]]))
                      sampVar <- factor(sapply(contigs[,sample], function(x){ x[1] }))
                  }else{
                      stopifnot(length(sample) == length(contigs))
                      sampVar <- factor(sample)
                  }
                  clono.list <- lapply(levels(sampVar), function(lv){
                      clono <- tenXquant(contigs[which(sampVar == lv)],
                                         sample = NULL)
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
              contigs <- contigs[contigs[,'productive']=='True']
              contigs <- contigs[contigs[,'chain'] %in% c('TRA','TRB')]
              
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
              clono <- clono[, colSums(clono) > 0]
              rownames(clono) <- names(contigs)
              
              return(clono)
          })

