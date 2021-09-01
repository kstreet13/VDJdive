


EMquant <- function(sce, TCRcol = 'contigs', thresh = .01, iter.max = 1000){
    contigs <- sce[[TCRcol]]
    
    # explore possible numbers of alpha and beta chains
    nAlpha <- sapply(contigs, function(x){
        sum(x$chain == 'TRA' & x$productive == 'True')
    })
    nBeta <- sapply(contigs, function(x){
        sum(x$chain == 'TRB' & x$productive == 'True')
    })
    
    # remove unproductive and 'Multi' contigs (for now)
    contigs <- contigs[contigs[,'productive']=='True']
    contigs <- contigs[contigs[,'chain'] %in% c('TRA','TRB')]
    
    # find all unique alpha chains
    all.alphas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRA']))
    
    # find all unique beta chains
    all.betas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']=='TRB']))
    
    # initialize counts matrix (#alpha-by-#beta)
    require(Matrix)
    counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), sparse = TRUE)
    rownames(counts) <- all.alphas
    colnames(counts) <- all.betas
    
    # step 1: assign cells with 1 alpha, 1 beta
    ############################################
    temp <- contigs[nAlpha == 1 & nBeta == 1]
    uniquecounts <- unclass(table(unlist(temp[temp[,'chain']=='TRA','cdr3']), 
                                  unlist(temp[temp[,'chain']=='TRB','cdr3'])))
    counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
    #rm(uniquecounts)
    
    # step 2: assign (multi-mapped) cells with >1 alpha or beta (equally across all 2-4 possibilities)
    ############################################################
    temp <- contigs[nAlpha > 0 & nBeta > 0 & nAlpha+nBeta > 2]
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
    temp <- contigs[(nAlpha == 0 | nBeta == 0) & nAlpha+nBeta > 0]
    for(i in seq_along(temp)){
        if(sum(temp[[i]]$chain=='TRA')==0){
            # only beta chain(s)
            rowind <- which(rowSums(counts.old[,temp[[i]]$cdr3,drop=FALSE]) > 0)
            if(length(rowind)==0){
                rowind <- seq_len(nrow(counts))
            }
            nposs <- nrow(temp[[i]]) * length(rowind)
            counts[rowind, temp[[i]]$cdr3] <- counts[rowind, temp[[i]]$cdr3] + 1/nposs
        }
        if(sum(temp[[i]]$chain=='TRB')==0){
            # only alpha chain(s)
            colind <- which(colSums(counts.old[temp[[i]]$cdr3,,drop=FALSE]) > 0)
            if(length(colind)==0){
                colind <- seq_len(ncol(counts))
            }
            nposs <- nrow(temp[[i]]) * length(colind)
            counts[temp[[i]]$cdr3, colind] <- counts[temp[[i]]$cdr3, colind] + 1/nposs
        }
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
        temp <- contigs[nAlpha > 0 & nBeta > 0 & nAlpha+nBeta > 2]
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
        temp <- contigs[(nAlpha == 0 | nBeta == 0) & nAlpha+nBeta > 0]
        for(i in seq_along(temp)){
            if(sum(temp[[i]]$chain=='TRA')==0){
                # only beta chain(s)
                rowind <- which(rowSums(counts.old[,temp[[i]]$cdr3,drop=FALSE]) > 0)
                replacement <- counts[rowind, temp[[i]]$cdr3] +
                    counts.old[rowind, temp[[i]]$cdr3] /
                    sum(counts.old[rowind, temp[[i]]$cdr3])
                counts[rowind, temp[[i]]$cdr3] <- replacement
            }
            if(sum(temp[[i]]$chain=='TRB')==0){
                # only alpha chain(s)
                colind <- which(colSums(counts.old[temp[[i]]$cdr3,,drop=FALSE]) > 0)
                replacement <- counts[temp[[i]]$cdr3, colind] +
                    counts.old[temp[[i]]$cdr3, colind] /
                    sum(counts.old[temp[[i]]$cdr3, colind])
                counts[temp[[i]]$cdr3, colind] <- replacement
            }
        }
        
        # check for convergence
        diff <- sum(abs(counts-counts.old))
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
    
    # update sce
    colData(sce)$clono <- clono
    return(sce)
    # just return vector of non-zero values, for diversity calculations
    #return(counts@x)
}

