# clonoStats helper functions

# remove all non-essential contig data
# (ie. for TCRs, only keep alpha and beta chains)
.prepContigs <- function(contigs, type){
    # remove unproductive and 'Multi' contigs (for now?)
    contigs <- contigs[contigs[,'productive']]
    if(type == 'TCR'){
        type1 <- 'TRA'
        type2 <- 'TRB'
    }
    if(type == 'BCR'){
        type1 <- 'IGH'
        type2 <- c('IGL','IGK')
        # prepend "IGK" or "IGL" to distinguish
        cellbarcodes <- names(contigs)
        contigs <- unlist(contigs)
        t2ind <- which(contigs$chain %in% type2)
        contigs$cdr3[t2ind] <- paste(contigs$chain[t2ind], 
                                     contigs$cdr3[t2ind], sep = '-')
        contigs <- split(DataFrame(contigs), 
                         factor(contigs$barcode, cellbarcodes))
    }
    contigs <- contigs[contigs[,'chain'] %in% c(type1, type2)]
    return(contigs)
}

######################################
# single-sample clonotype assignment #
######################################
# EM algorithm
.EM_sample <- function(contigs, type, lang, thresh, iter.max){
    contigs <- .prepContigs(contigs, type)
    
    if(type == 'TCR'){
        type1 <- 'TRA'
        type2 <- 'TRB'
    }
    if(type == 'BCR'){
        type1 <- 'IGH'
        type2 <- c('IGL','IGK')
    }
    
    # numbers of alpha and beta chains per cell
    nAlpha <- sum(contigs[,"chain"] %in% type1)
    nBeta <- sum(contigs[,"chain"] %in% type2)
    
    # find all unique alpha chains
    all.alphas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']%in%type1]))
    all.alphas <- c('UNKNOWN_1', all.alphas)
    # find all unique beta chains
    all.betas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']%in%type2]))
    all.betas <- c('UNKNOWN_2', all.betas)

    # initialize counts matrix (#alpha-by-#beta)
    counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), 
                     sparse = TRUE)
    rownames(counts) <- all.alphas
    colnames(counts) <- all.betas
    
    # useful indices
    ind.unique <- which(nAlpha == 1 & nBeta == 1)
    ind.multiple <- which(nAlpha > 0 & nBeta > 0 & nAlpha+nBeta > 2)
    ind.noAlpha <- which(nAlpha == 0 & nBeta > 0)
    ind.noBeta <- which(nBeta == 0 & nAlpha > 0)
    ind.ambiguous <- sort(c(ind.multiple, ind.noAlpha, ind.noBeta))
    
    # use which alpha sequence and which beta sequence each contig represents to
    # calculate its (possible) index in the clonotype matrix
    wa <- match(contigs[,'cdr3'][contigs[,'chain']%in%type1], all.alphas)
    wb <- match(contigs[,'cdr3'][contigs[,'chain']%in%type2], all.betas)
    # this takes care of the 1-alpha + 1-beta indices without looping
    poss.indices <- as.list(length(all.alphas)*(wb-1L) + wa)
    # others are incorrect, though
    poss.indices[ind.multiple] <- lapply(ind.multiple, function(i){
        return(sort(unique(as.integer(outer(wa[[i]], wb[[i]], 
                                            FUN = function(a,b){
                                                length(all.alphas)*(b-1) + a
                                            })))))
    })
    poss.indices[ind.noAlpha] <- lapply(ind.noAlpha, function(i){
        a.i <- seq_along(all.alphas)
        return(sort(unique(as.integer(outer(a.i, wb[[i]], 
                                            FUN = function(a,b){
                                                length(all.alphas)*(b-1) + a
                                            })))))
    })
    poss.indices[ind.noBeta] <- lapply(ind.noBeta, function(i){
        b.i <- seq_along(all.betas)
        return(sort(unique(as.integer(outer(wa[[i]], b.i, 
                                            FUN = function(a,b){
                                                length(all.alphas)*(b-1) + a
                                            })))))
    })
    
    # step 1: assign cells with 1 alpha, 1 beta
    ############################################
    temp <- contigs[ind.unique]
    uniquecounts <- unclass(table(unlist(temp[temp[,'chain']%in%type1,'cdr3']),
                                  unlist(temp[temp[,'chain']%in%type2,'cdr3'])))
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
    if(lang == 'python'){
        # setup inputs
        t.indices <- poss.indices[ind.ambiguous]
        # so python takes it as a list of lists, not a dict:
        names(t.indices) <- NULL
        # force python to take length-1 vectors as lists
        l1.idx <- which(lengths(t.indices) == 1)
        for(ii in l1.idx){
            t.indices[[ii]] <- list(as.integer(t.indices[[ii]]))
        }
        
        # iteration handled by python, via basilisk
        cl <- basiliskStart(pyenv)
        counts <- basiliskRun(cl, function(uniquecounts, counts.old, t.indices, 
                                           thresh, iter.max){
            mod <- reticulate::import(module = "vdjHelpers", 
                                      convert = TRUE)
            return(mod$TCR_EM_counts(uniquecounts, counts.old, t.indices, 
                                     thresh, iter.max))
        }, uniquecounts = uniquecounts, counts.old = counts.old, 
        t.indices = t.indices, thresh = thresh, iter.max = iter.max)
        basiliskStop(cl)
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
    clonoP <- as.integer(c(0,cumsum(lengths(poss.indices))))
    clonoX <- unlist(lapply(which(lengths(poss.indices) > 0), function(i){
        counts[poss.indices[[i]]] /
            sum(counts[poss.indices[[i]]])
    }))
    clono <- new('dgRMatrix', j = clonoJ, p = clonoP, x = clonoX,
                 Dim = as.integer(c(length(contigs), length(counts))))
    colnames(clono) <- as.character(outer(all.alphas, all.betas,
                                          FUN = paste))
    clono <- clono[, which(colSums(clono) > 0), drop = FALSE]
    rownames(clono) <- names(contigs)
    
    return(clono)
}

# CellRanger labels
.CR_sample <- function(contigs, type){
    contigs <- .prepContigs(contigs, type)    
    
    clonoID <- vapply(contigs, function(x){
        alphas <- sort(x$cdr3[x$chain=='TRA'])
        betas <- sort(x$cdr3[x$chain=='TRB'])
        return(paste(paste(alphas, collapse = ','),
                     paste(betas, collapse = ',')))
    }, FUN.VALUE = 'A')
    
    all_clonotypes <- unique(clonoID)
    
    poss.indices <- match(clonoID, all_clonotypes)
    
    
    # build full cells-by-clonotypes matrix
    clonoJ <- as.integer(unlist(poss.indices)-1)
    clonoP <- as.integer(c(0,cumsum(lengths(poss.indices))))
    clonoX <- rep(1, sum(lengths(poss.indices)))
    clono <- new('dgRMatrix', j = clonoJ, p = clonoP, x = clonoX,
                 Dim = as.integer(c(length(contigs), length(all_clonotypes))))
    colnames(clono) <- all_clonotypes
    clono <- clono[, which(colSums(clono) > 0), drop = FALSE]
    rownames(clono) <- names(contigs)
    # remove empty clonotype (called "None" in some versions)
    clono <- clono[,-which(colnames(clono) == " ")]
    
    return(clono)
}

# Unique assignment
.UNIQ_sample <- function(contigs, type){
    contigs <- .prepContigs(contigs, type)
    
    if(type == 'TCR'){
        type1 <- 'TRA'
        type2 <- 'TRB'
    }
    if(type == 'BCR'){
        type1 <- 'IGH'
        type2 <- c('IGL','IGK')
    }
    
    # numbers of alpha and beta chains per cell
    nAlpha <- sum(contigs[,"chain"] %in% type1)
    nBeta <- sum(contigs[,"chain"] %in% type2)
    
    # find all unique alpha chains
    all.alphas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']%in%type1]))
    if(length(all.alphas)==0){
        # can't assign clonotypes
        # return empty matrix, in case this is just one sample of many
        clono <- Matrix(0, nrow = length(contigs), ncol = 0)
        rownames(clono) <- names(contigs)
        return(clono)
    }
    # find all unique beta chains
    all.betas <- unique(unlist(contigs[,'cdr3'][contigs[,'chain']%in%type2]))
    if(length(all.betas)==0){
        # can't assign clonotypes
        # return empty matrix, in case this is just one sample of many
        clono <- Matrix(0, nrow = length(contigs), ncol = 0)
        rownames(clono) <- names(contigs)
        return(clono)
    }
    
    # initialize counts matrix (#alpha-by-#beta)
    counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), 
                     sparse = TRUE)
    rownames(counts) <- all.alphas
    colnames(counts) <- all.betas
    
    # useful indices
    ind.unique <- which(nAlpha == 1 & nBeta == 1)
    
    # use which alpha sequence and which beta sequence each contig represents 
    # to calculate its (possible) index in the clonotype matrix
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
}


.CHN_sample <- function(contigs, method){
    # in this case, 'method' indicates the type of chain being quantified
    contigs <- contigs[contigs[,'productive']]
    contigs <- contigs[contigs[,'chain'] == method]
    
    # find all unique chains (not necessarily alphas)
    all.alphas <- unique(unlist(contigs[,'cdr3']))
    if(length(all.alphas)==0){
        # can't assign clonotypes
        # return empty matrix, in case this is just one sample of many
        clono <- Matrix(0, nrow = length(contigs), ncol = 0)
        rownames(clono) <- names(contigs)
        return(clono)
    }
    
    # build full cells-by-clonotypes matrix
    clono <- Matrix(0, nrow = length(contigs), ncol = length(all.alphas))
    for(i in seq_along(contigs)){
        clono[i,] <- as.numeric(all.alphas %in% contigs[[i]][,'cdr3'])
    }
    colnames(clono) <- all.alphas
    rownames(clono) <- names(contigs)
    
    return(clono)
}




