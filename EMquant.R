require(SingleCellExperiment)
sce <- readRDS('~/Desktop/toyTCRdata.rds')

# explore possible numbers of alpha and beta chains
nAlpha <- sapply(sce$contigs, function(x){
    sum(x$chain == 'TRA' & x$productive == 'True')
})
nBeta <- sapply(sce$contigs, function(x){
    sum(x$chain == 'TRB' & x$productive == 'True')
})
table(nAlpha,nBeta)
# 

# remind myself how this works
sce$contigs[[4]]
sce$contigs[,'umis']
head(unlist(sce$contigs[,'umis']))
is.prod <- sce$contigs[,"productive"]=="True"
has.prod <- any(is.prod)


# remove unproductive and 'Multi' contigs (for now)
sce$contigs <- sce$contigs[sce$contigs[,'productive']=='True']
sce$contigs <- sce$contigs[sce$contigs[,'chain'] %in% c('TRA','TRB')]

# find all unique alpha chains
all.alphas <- unique(unlist(sce$contigs[,'cdr3'][sce$contigs[,'chain']=='TRA']))

# find all unique beta chains
all.betas <- unique(unlist(sce$contigs[,'cdr3'][sce$contigs[,'chain']=='TRB']))

# initialize counts matrix
require(Matrix)
counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), sparse = TRUE)
rownames(counts) <- all.alphas
colnames(counts) <- all.betas


# step 1: assign cells with 1 alpha, 1 beta
############################################
temp <- sce$contigs[nAlpha == 1 & nBeta == 1]
uniquecounts <- unclass(table(unlist(temp[temp[,'chain']=='TRA','cdr3']), 
                              unlist(temp[temp[,'chain']=='TRB','cdr3'])))
counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
#rm(uniquecounts)

# step 2: assign (multi-mapped) cells with >1 alpha or beta (proportionally across all 2-4 possibilities)
############################################################
temp <- sce$contigs[nAlpha > 0 & nBeta > 0 & nAlpha+nBeta > 2]
for(i in seq_along(temp)){
    nposs <- sum(temp[[i]]$chain=='TRA') * sum(temp[[i]]$chain=='TRB')
    replacement <- counts[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
                          temp[[i]]$cdr3[temp[[i]]$chain=='TRB']] + 1/nposs
    counts[temp[[i]]$cdr3[temp[[i]]$chain=='TRA'],
           temp[[i]]$cdr3[temp[[i]]$chain=='TRB']] <- replacement
}
counts.old <- counts

# step 3: assign (multi-mapped) cells with 0 alpha or beta (proportionally across all 1-? possibilities)
###########################################################
temp <- sce$contigs[(nAlpha == 0 | nBeta == 0) & nAlpha+nBeta > 0]
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
diff <- NULL
while(working){
    counts <- Matrix(0, nrow = length(all.alphas), ncol = length(all.betas), sparse = TRUE)
    rownames(counts) <- all.alphas
    colnames(counts) <- all.betas
    counts[rownames(uniquecounts), colnames(uniquecounts)] <- uniquecounts
    
    # step 2
    #########
    temp <- sce$contigs[nAlpha > 0 & nBeta > 0 & nAlpha+nBeta > 2]
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
    
    # step 3
    #########
    temp <- sce$contigs[(nAlpha == 0 | nBeta == 0) & nAlpha+nBeta > 0]
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
            nposs <- nrow(temp[[i]]) * length(colind)
            replacement <- counts[temp[[i]]$cdr3, colind] +
                counts.old[temp[[i]]$cdr3, colind] /
                sum(counts.old[temp[[i]]$cdr3, colind])
            counts[temp[[i]]$cdr3, colind] <- replacement
        }
    }
    
    diff <- c(diff,sum(abs(counts-counts.old)))
    print(diff[length(diff)])
    if(diff[length(diff)] < .1){
        working <- FALSE
    }
    
    counts.old <- counts
}
