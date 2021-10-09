

#' @param sce SingleCellExperiment object with TCR data
#' @param TCRcol element of the colData of sce that contains the TCR data
#' @param sample optional variable for splitting quantification by sample. This can speed up computation and avoid crosstalk between samples.
#' @param thresh threshold for convergence in assigning clonotype counts
#' @param iter.max maximum number of iterations for the EM algorithm.


EMquant <- function(sce, TCRcol = 'contigs', sample = NULL, thresh = .01, iter.max = 1000){
    if(!is.null(sample)){
        if(length(sample) == 1){
            stopifnot(sample %in% names(colData(sce)))
            sampVar <- factor(sce[[sample]])
        }else{
            stopifnot(length(sample) == ncol(sce))
            sampVar <- factor(sample)
        }
        clono.list <- lapply(levels(sampVar), function(lv){
            EMquant(sce[,which(sampVar == lv)],
                    TCRcol = TCRcol, sample = NULL,
                    thresh = thresh, iter.max = iter.max)$clono
        })
        all.clonotypes <- unique(unlist(lapply(clono.list, colnames)))
        clono.list <- lapply(clono.list, function(mat){
            missing <- all.clonotypes[! all.clonotypes %in% colnames(mat)]
            zeros <- Matrix(0, ncol = length(missing), nrow = nrow(mat), 
                            sparse = TRUE)
            colnames(zeros) <- missing
            return(cbind(mat, zeros))
        })
        clono <- do.call(rbind, clono.list)
        clono <- clono[colnames(sce), ]
        # update sce
        colData(sce)$clono <- clono
        return(sce)
    }
    
    contigs <- sce[[TCRcol]]
    
    # remove unproductive and 'Multi' contigs (for now)
    contigs <- contigs[contigs[,'productive']=='True']
    contigs <- contigs[contigs[,'chain'] %in% c('TRA','TRB')]
    
    # explore possible numbers of alpha and beta chains
    nAlpha <- sum(contigs[,"chain"] == 'TRA')
    nBeta <- sum(contigs[,"chain"] == 'TRB')
    
    # find all unique alpha chains
    all.alphas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRA']))
    
    # find all unique beta chains
    all.betas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRB']))
    
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
    
    # step 1: assign cells with 1 alpha, 1 beta
    ############################################
    temp <- contigs[ind.unique]
    uniquecounts <- unclass(table(unlist(temp[temp[,'chain']=='TRA','cdr3']), 
                                  unlist(temp[temp[,'chain']=='TRB','cdr3'])))
    counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
    #rm(uniquecounts)
    
    # step 2: assign (multi-mapped) cells with >1 alpha or beta (equally across all 2-4 possibilities)
    ############################################################
    temp <- contigs[ind.multiple]
    for(i in seq_along(temp)){
        nposs <- sum(temp[[i]]$chain=='TRA') * sum(temp[[i]]$chain=='TRB')
        replacement <- counts[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
                              temp[[i]]$cdr3[temp[[i]]$chain=='TRB']] + 1/nposs
        counts[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
               temp[[i]]$cdr3[temp[[i]]$chain=='TRB']] <- replacement
    }
    counts.old <- counts
    
    # step 3: assign (multi-mapped) cells with 0 alpha or beta (equally across all 1-? possibilities)
    ###########################################################
    temp <- contigs[ind.noAlpha]
    for(i in seq_along(temp)){
        # only beta chain(s)
        rowind <- which(rowSums(counts.old[,temp[[i]]$cdr3,drop=FALSE]) > 0)
        if(length(rowind)==0){
            rowind <- seq_len(nrow(counts))
        }
        nposs <- nrow(temp[[i]]) * length(rowind)
        counts[rowind, temp[[i]]$cdr3] <- counts[rowind, temp[[i]]$cdr3] + 1/nposs
    }
    temp <- contigs[ind.noBeta]
    for(i in seq_along(temp)){
        # only alpha chain(s)
        colind <- which(colSums(counts.old[temp[[i]]$cdr3,,drop=FALSE]) > 0)
        if(length(colind)==0){
            colind <- seq_len(ncol(counts))
        }
        nposs <- nrow(temp[[i]]) * length(colind)
        counts[temp[[i]]$cdr3, colind] <- counts[temp[[i]]$cdr3, colind] + 1/nposs
    }
    counts.old <- counts
    
    # repeat 2 and 3 (proportional to previous counts)
    #################
    working <- TRUE
    iters <- 0
    while(working){
        iters <- iters + 1
        counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), sparse = TRUE)
        rownames(counts) <- all.alphas
        colnames(counts) <- all.betas
        counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
        
        # step 2 (proportional to counts.old)
        #########
        temp <- contigs[ind.multiple]
        for(i in seq_along(temp)){
            replacement <- counts[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
                                  temp[[i]]$cdr3[temp[[i]]$chain=='TRB']] + 
                counts.old[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
                           temp[[i]]$cdr3[temp[[i]]$chain=='TRB']] / 
                sum(counts.old[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
                               temp[[i]]$cdr3[temp[[i]]$chain=='TRB']])
            counts[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
                   temp[[i]]$cdr3[temp[[i]]$chain=='TRB']] <- replacement
        }
        
        # step 3 (proportional to counts.old)
        #########
        temp <- contigs[ind.noAlpha]
        for(i in seq_along(temp)){
            # only beta chain(s)
            rowind <- which(rowSums(counts.old[,temp[[i]]$cdr3,drop=FALSE]) > 0)
            counts[rowind, temp[[i]]$cdr3] <- 
                counts[rowind, temp[[i]]$cdr3] +
                counts.old[rowind, temp[[i]]$cdr3] /
                sum(counts.old[rowind, temp[[i]]$cdr3])
        }
        temp <- contigs[ind.noBeta]
        for(i in seq_along(temp)){
            # only alpha chain(s)
            colind <- which(colSums(counts.old[temp[[i]]$cdr3,,drop=FALSE]) > 0)
            counts[temp[[i]]$cdr3, colind] <- 
                counts[temp[[i]]$cdr3, colind] +
                counts.old[temp[[i]]$cdr3, colind] /
                sum(counts.old[temp[[i]]$cdr3, colind])
        }
        
        # check for convergence
        diff <- max(abs(counts-counts.old))
        #print(diff[length(diff)])
        if(diff < thresh | iters >= iter.max){
            working <- FALSE
            # build full cells-by-clonotypes matrix
            clono <- Matrix(0, nrow = ncol(sce), ncol = length(all.alphas)*length(all.betas), sparse = TRUE)
            colnames(clono) <- as.character(outer(rownames(counts), 
                                                  colnames(counts), 
                                                  FUN = paste))
            for(i in seq_along(contigs)){
                if(nAlpha[i] == 0 | nBeta[i] == 0){
                    if(nAlpha[i] + nBeta[i] == 0){
                        next
                    }else{
                        if(nAlpha[i] == 0){
                            # only beta chain(s)
                            rowind <- which(rowSums(counts.old[,contigs[[i]]$cdr3,drop=FALSE]) > 0)
                            contribution <- counts.old[rowind, contigs[[i]]$cdr3, drop = FALSE] /
                                sum(counts.old[rowind, contigs[[i]]$cdr3])
                        }
                        if(nBeta[i] == 0){
                            # only alpha chain(s)
                            colind <- which(colSums(counts.old[contigs[[i]]$cdr3,,drop=FALSE]) > 0)
                            contribution <- counts.old[contigs[[i]]$cdr3, colind, drop = FALSE] /
                                sum(counts.old[contigs[[i]]$cdr3, colind])
                        }
                    }
                }else{
                    contribution <- counts.old[contigs[[i]]$cdr3[contigs[[i]]$chain=='TRA'],
                                               contigs[[i]]$cdr3[contigs[[i]]$chain=='TRB'],
                                               drop = FALSE] / 
                        sum(counts.old[contigs[[i]]$cdr3[contigs[[i]]$chain=='TRA'],
                                       contigs[[i]]$cdr3[contigs[[i]]$chain=='TRB']])
                }
                if(length(contribution) > 0){
                    col.ind <- outer(match(rownames(contribution), rownames(counts)),
                                     (match(colnames(contribution), colnames(counts))-1) * nrow(counts),
                                     FUN = "+")
                    clono[i, as.numeric(col.ind)] <- as.numeric(contribution)
                }
            }
            clono <- clono[, colSums(clono) > 0]
        }
        
        counts.old <- counts
    }
    rownames(clono) <- colnames(sce)
    
    # update sce
    colData(sce)$clono <- clono
    return(sce)
    # just return vector of non-zero values, for diversity calculations
    #return(counts@x)
}




# this is the fast version
# results are slightly different, need to figure out why
EMquant2 <- function(sce, TCRcol = 'contigs', sample = NULL, thresh = .01, iter.max = 1000){
    if(!is.null(sample)){
        if(length(sample) == 1){
            stopifnot(sample %in% names(colData(sce)))
            sampVar <- factor(sce[[sample]])
        }else{
            stopifnot(length(sample) == ncol(sce))
            sampVar <- factor(sample)
        }
        clono.list <- lapply(levels(sampVar), function(lv){
            EMquant2(sce[,which(sampVar == lv)],
                     TCRcol = TCRcol, sample = NULL,
                     thresh = thresh, iter.max = iter.max)$clono
        })
        all.clonotypes <- unique(unlist(lapply(clono.list, colnames)))
        clono.list <- lapply(clono.list, function(mat){
            missing <- all.clonotypes[! all.clonotypes %in% colnames(mat)]
            zeros <- Matrix(0, ncol = length(missing), nrow = nrow(mat), 
                            sparse = TRUE)
            colnames(zeros) <- missing
            return(cbind(mat, zeros))
        })
        clono <- do.call(rbind, clono.list)
        clono <- clono[colnames(sce), ]
        # update sce
        colData(sce)$clono <- clono
        return(sce)
    }
    
    contigs <- sce[[TCRcol]]
    
    # remove unproductive and 'Multi' contigs (for now?)
    contigs <- contigs[contigs[,'productive']=='True']
    contigs <- contigs[contigs[,'chain'] %in% c('TRA','TRB')]
    
    # explore possible numbers of alpha and beta chains
    nAlpha <- sum(contigs[,"chain"] == 'TRA')
    nBeta <- sum(contigs[,"chain"] == 'TRB')
    
    # find all unique alpha chains
    all.alphas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRA']))
    
    # find all unique beta chains
    all.betas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRB']))
    
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
    poss.indices <- as.list(length(all.alphas)*(wb-1) + wa)
    # others are incorrect, though
    poss.indices[ind.multiple] <- sapply(ind.multiple, function(i){
        return(as.numeric(sapply(wa[[i]], function(a){
            sapply(wb[[i]], function(b){
                length(all.alphas)*(b-1) + a
            })
        })))
    })
    poss.indices[ind.noAlpha] <- sapply(ind.noAlpha, function(i){
        a.i <- seq_along(all.alphas)
        return(as.numeric(sapply(wb[[i]], function(b){
            length(all.alphas)*(b-1) + a.i
        })))
    })
    poss.indices[ind.noBeta] <- sapply(ind.noBeta, function(i){
        b.i <- seq_along(all.betas)
        return(as.numeric(sapply(wa[[i]], function(a){
            length(all.alphas)*(b.i-1) + a
        })))
    })
    
    # step 1: assign cells with 1 alpha, 1 beta
    ############################################
    temp <- contigs[ind.unique]
    uniquecounts <- unclass(table(unlist(temp[temp[,'chain']=='TRA','cdr3']), 
                                  unlist(temp[temp[,'chain']=='TRB','cdr3'])))
    counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
    uniquecounts <- counts <- as.numeric(counts)
    
    # step 2: assign (multi-mapped) cells (equally across all possibilities)
    ############################################################
    temp <- contigs[ind.ambiguous]
    t.indices <- poss.indices[ind.ambiguous]
    for(i in seq_along(temp)){
        counts[t.indices[[i]]] <- 
            counts[t.indices[[i]]] + 1/length(t.indices[[i]])
    }
    counts.old <- counts
    
    # repeat 2 (proportional to previous counts)
    #################
    working <- TRUE
    iters <- 0
    while(working){
        iters <- iters + 1
        counts <- uniquecounts
        
        # step 2 (proportional to counts.old)
        #########
        temp <- contigs[ind.ambiguous]
        t.indices <- poss.indices[ind.ambiguous]
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
            # build full cells-by-clonotypes matrix
            clono <- Matrix(0, nrow = ncol(sce), 
                            ncol = length(all.alphas)*length(all.betas), 
                            sparse = TRUE)
            colnames(clono) <- as.character(outer(all.alphas, all.betas, 
                                                  FUN = paste))
            
            for(i in ind.unique){
                clono[i, poss.indices[[i]]] <- 1
            }
            for(i in ind.ambiguous){
                # use counts.old so that it's equal to the current contribution, not running one extra iteration
                clono[i, poss.indices[[i]]] <- 
                    counts.old[poss.indices[[i]]] / 
                    sum(counts.old[poss.indices[[i]]])
            }
            
            clono <- clono[, colSums(clono) > 0]
        }
        
        counts.old <- counts
    }
    rownames(clono) <- colnames(sce)
    
    # update sce
    colData(sce)$clono <- clono
    return(sce)
    # just return vector of non-zero values, for diversity calculations
    #return(counts@x)
}

