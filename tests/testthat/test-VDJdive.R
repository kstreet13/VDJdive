context("Test VDJdive.")
library(utils) # needed for data()
library(stats) # for rpois()
library(S4Vectors) # for metadata()

test_that("utility functions work", {
    # load example data
    data("contigs")
    
    # make SCE object with matching barcodes and sample IDs
    ncells <- 24
    u <- matrix(rpois(1000 * ncells, 5), ncol = ncells)
    barcodes <- vapply(contigs[,'barcode'], function(x){ x[1] }, 'A')
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = u),
        colData = data.frame(Barcode = barcodes,
                             group = samples))
    sce <- addVDJtoSCE(contigs, sce)
    
    # errors before running clonoStats
    expect_error(splitClonotypes(sce, by = samples), 
                 'No clonotype counts found')
    expect_error(summarizeClonotypes(sce, by = samples), 
                 'No clonotype counts found')
    
    sce <- clonoStats(sce, method = 'unique', assignment = TRUE)
    
    # splitClonotypes (CL = counts list)
    CL <- splitClonotypes(sce, by = 'sample')
    expect_equal(length(CL), 2)
    expect_equivalent(dim(CL[[1]]), c(12, 7))
    expect_equivalent(dim(CL[[2]]), c(12, 7))
    expect_equal(sum(CL[[1]]), 6)
    
    CL2 <- splitClonotypes(as.matrix(metadata(sce)$clonoStats@assignment), 
                           by = samples)
    expect_identical(CL, CL2)
    
    
    # summarizeClonotypes (SLC = sample-level counts)
    SLC <- summarizeClonotypes(sce, by = 'sample')
    expect_equivalent(dim(SLC), c(7, 2))
    expect_equivalent(colSums(SLC), c(6, 6))
    
    SLC2 <- summarizeClonotypes(as.matrix(metadata(sce)$clonoStats@assignment),
                                by = samples)
    expect_identical(SLC, SLC2)
    
    SLF <- summarizeClonotypes(sce, by = 'sample', mode = 'tab')
    expect_equivalent(dim(SLF), c(4, 2))
    expect_equivalent(rowSums(SLF), c(5,7,1,1))
    
    cn <- clonoNames(sce)
    expect_equal(length(unlist(strsplit(cn,split=' '))), 2 * length(cn))
})

test_that("input/output functions work", {
    # load example data
    data("contigs")

    # write the example data to a temporary directory
    loc <- tempdir()
    writeVDJcontigs(loc, contigs)
    expect_true(file.exists(file.path(loc,'sample1',
                                      'filtered_contig_annotations.csv')))

    # specify sample locations and read in data
    samples <- file.path(loc, c('sample1','sample2'))
    contigs <- readVDJcontigs(samples)
    expect_equivalent(lengths(contigs),
                      c(3,6,2,4,2,3,2,2,2,4,2,2,3,1,2,4,2,2,1,2,2,1,1,3))

    # make SCE object with matching barcodes and sample IDs
    ncells <- 24
    u <- matrix(rpois(1000 * ncells, 5), ncol = ncells)
    barcodes <- vapply(contigs[,'barcode'], function(x){ x[1] }, 'A')
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = u),
        colData = data.frame(Barcode = barcodes,
                             group = samples))

    sce <- addVDJtoSCE(contigs, sce)
    expect_true('contigs' %in% names(colData(sce)))
    expect_equivalent(lengths(sce$contigs),
                      c(3,6,2,4,2,3,2,2,2,4,2,2,3,1,2,4,2,2,1,2,2,1,1,3))
    
    # from file, with some loss
    sce$contigs <- NULL
    samples <- file.path(loc, c('sample1','sample2'))
    sce <- sce[, -1]
    expect_message({sce <- addVDJtoSCE(samples, sce)}, 
                   regexp = '1 cells with V')
    expect_equivalent(lengths(sce$contigs),
                      c(6,2,4,2,3,2,2,2,4,2,2,3,1,2,4,2,2,1,2,2,1,1,3))
})

test_that("clonoStats function works as expected", {
    # load example data
    data("contigs")

    uniq <- clonoStats(contigs, method = 'unique', assignment = TRUE)
    expect_is(uniq, 'clonoStats')
    expect_equivalent(dim(uniq@abundance), c(7,2))
    expect_equal(sum(uniq@abundance), 12)
    expect_equivalent(dim(uniq@frequency), c(4,2))
    expect_equal(sum(uniq@frequency), 14)
    expect_equivalent(dim(uniq@assignment), c(24, 7))
    expect_equal(sum(uniq@assignment), 12)

    crng <- clonoStats(contigs, method = 'CellRanger', assignment = TRUE)
    expect_is(crng, 'clonoStats')
    expect_equivalent(dim(crng@abundance), c(17,2))
    expect_equal(sum(crng@abundance), 22)
    expect_equivalent(dim(crng@frequency), c(4,2))
    expect_equal(sum(crng@frequency), 34)
    expect_equivalent(dim(crng@assignment), c(24, 17))
    expect_equal(sum(crng@assignment), 22)

    alph <- clonoStats(contigs, method = 'TRA', assignment = TRUE)
    expect_is(alph, 'clonoStats')
    expect_equivalent(dim(alph@abundance), c(18,2))
    expect_equal(sum(alph@abundance), 24)
    expect_equivalent(dim(alph@frequency), c(4,2))
    expect_equal(sum(alph@frequency), 36)
    expect_equivalent(dim(alph@assignment), c(24, 18))
    expect_equal(sum(alph@assignment), 24)
    expect_equivalent(grep('\\s$', clonoNames(alph)), seq_len(18))
    
    emal <- clonoStats(contigs, method = 'EM', assignment = TRUE)
    expect_is(emal, 'clonoStats')
    expect_equivalent(dim(emal@abundance), c(60,2))
    expect_equal(sum(emal@abundance), 22)
    expect_equivalent(dim(emal@frequency), c(4,2))
    expect_equal(sum(emal@frequency), 120)
    expect_equivalent(dim(emal@assignment), c(24, 60))
    expect_equal(sum(emal@assignment), 22)
    
    emal2 <- clonoStats(contigs, method = 'EM', 
                        assignment = TRUE, lang = 'r')
    expect_is(emal2, 'clonoStats')
    expect_equivalent(dim(emal2@abundance), c(60,2))
    expect_equal(sum(emal2@abundance), 22)
    expect_equivalent(dim(emal2@frequency), c(4,2))
    expect_equal(sum(emal2@frequency), 120)
    expect_equivalent(dim(emal2@assignment), c(24, 60))
    expect_equal(sum(emal2@assignment), 22)
    expect_true(max(abs(emal2@assignment - emal@assignment)) < .0001)

    # make SCE object with matching barcodes and sample IDs
    ncells <- 24
    u <- matrix(rpois(1000 * ncells, 5), ncol = ncells)
    barcodes <- vapply(contigs[,'barcode'], function(x){ x[1] }, 'A')
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = u),
        colData = data.frame(Barcode = barcodes,
                             group = samples))
    sce <- addVDJtoSCE(contigs, sce)
    
    sceUniq <- clonoStats(sce, method = 'unique', assignment = TRUE)
    expect_identical(uniq, metadata(sceUniq)$clonoStats)
    
    sceCR <- clonoStats(sce, method = 'CellRanger', assignment = TRUE)
    expect_identical(crng, metadata(sceCR)$clonoStats)
    
    sceEM <- clonoStats(sce, method = 'EM', assignment = TRUE)
    expect_identical(emal, metadata(sceEM)$clonoStats)
    
    sceEMsamp <- clonoStats(sce, group = 'sample')
    sceEMsamp2 <- clonoStats(sce, group = sce$sample)
    expect_identical(metadata(sceEMsamp)$clonoStats, 
                     metadata(sceEMsamp2)$clonoStats)
    
    # clonoStats on object of type clonoStats
    clus <- sample(1:2, 24, replace = TRUE)
    emal3 <- clonoStats(emal, group = clus)
    expect_is(emal3, 'clonoStats')
    expect_equivalent(dim(emal3@abundance), c(60,2))
    expect_equal(sum(emal3@abundance), 22)
    expect_equal(ncol(emal3@frequency), 2)
    expect_equal(sum(emal3@frequency), 120)
    expect_equivalent(dim(emal3@assignment), c(24, 60))
    expect_equal(sum(emal3@assignment), 22)
    expect_true(max(abs(emal3@assignment - emal@assignment)) < .0001)
})

test_that("assignment functions handle edge cases", {
    # load example data
    data("contigs")
    
    # empty sample
    sample <- factor(vapply(contigs[,'sample'], function(x){ x[1] }, 'A'))
    levels(sample) <- c('sample1','sample2','sample3')
    uniq <- clonoStats(contigs, group = sample, method = 'unique')
    expect_is(uniq, 'clonoStats')
    expect_equivalent(dim(uniq@abundance), c(7,3))
    expect_equal(sum(uniq@abundance), 12)
    expect_equivalent(dim(uniq@frequency), c(4,3))
    expect_equal(sum(uniq@frequency), 21)
    
    # sample NAs
    sampNA <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    sampNA <- factor(sampNA, levels = c('sample1','sample3'))
    uniq <- clonoStats(contigs, group = sampNA, method = 'unique')
    expect_is(uniq, 'clonoStats')
    expect_equivalent(dim(uniq@abundance), c(3,2))
    expect_equal(sum(uniq@abundance), 6)
    expect_equivalent(dim(uniq@frequency), c(4,2))
    expect_equal(sum(uniq@frequency), 6)

    # SCE object with EXTRA CELLS and SAMPLE
    ncells <- 30
    u <- matrix(rpois(1000 * ncells, 5), ncol = ncells)
    barcodes <- vapply(contigs[,'barcode'], function(x){ x[1] }, 'A')
    barcodes <- c(barcodes, LETTERS[1:6])
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    samples <- c(samples, rep(c('sample1','sample2','sample3'), each = 2))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = u),
        colData = data.frame(Barcode = barcodes,
                             group = samples))
    sce <- addVDJtoSCE(contigs, sce)
    sceUniq <- clonoStats(sce, group ='sample', 
                          method = 'unique', assignment = TRUE)
    
    expect_is(metadata(sceUniq)$clonoStats, 'clonoStats')
    expect_equivalent(dim(metadata(sceUniq)$clonoStats@assignment), c(30, 7))
    expect_equal(sum(metadata(sceUniq)$clonoStats@assignment[
        sce$sample=='sample3',]), 0)
    
    # BCR data and mixtures
    uniq <- clonoStats(contigs, type = 'BCR', 
                       method = 'unique', assignment = TRUE)
    expect_equivalent(dim(uniq@assignment), c(24, 0))
    contigs[contigs[,'chain']=='TRA','chain'] <- 'IGH'
    uniq <- clonoStats(contigs, type = 'BCR', 
                       method = 'unique', assignment = TRUE)
    expect_equivalent(dim(uniq@assignment), c(24, 0))
    contigs[contigs[,'chain']=='TRB','chain'] <- 'IGL'
    
    # let it detect BCR
    uniq <- clonoStats(contigs, method = 'unique', assignment = TRUE)
    expect_is(uniq, 'clonoStats')
    expect_equivalent(dim(uniq@abundance), c(7,2))
    expect_equal(sum(uniq@abundance), 12)
    expect_equivalent(dim(uniq@frequency), c(4,2))
    expect_equal(sum(uniq@frequency), 14)
    expect_equivalent(dim(uniq@assignment), c(24, 7))
    expect_equal(sum(uniq@assignment), 12)

})

test_that("diversity calculation works", {
    # load example data
    data("contigs")
    
    xCR <- clonoStats(contigs, method = 'CellRanger')
    xEM <- clonoStats(contigs, method = 'EM')
    
    expect_warning({
        divCR <- calculateDiversity(xCR)
    }, regexp = 'use t=3')
    expect_equivalent(dim(divCR), c(10, 2))
    expect_equivalent(colnames(divCR), c('sample1', 'sample2'))

    expect_warning({
        divEM <- calculateDiversity(xEM)
    }, regexp = 'not valid with non-integer')
    expect_equivalent(dim(divEM), c(10, 2))
    expect_equivalent(colnames(divEM), c('sample1', 'sample2'))
})

test_that("PCA function works", {
    # load example data
    data("contigs")
    
    clono <- clonoStats(contigs)
    
    pca1 <- runVDJPCA(clono)
    expect_equivalent(dim(pca1$x), c(2, 2))
    
    pca2 <- runVDJPCA(clono, unit = 'clonotypes')
    expect_equivalent(dim(pca2$x), c(60, 2))
})

test_that("plotting functions work", {
    # load example data
    data("contigs")
    
    x <- clonoStats(contigs)
    p1 <- barVDJ(x)
    expect_equal(class(p1$layers[[1]]$geom)[1], 'GeomCol')
    
    p2 <- barVDJ(x, bySample = FALSE, 
                 title = 'bar plot', legend = TRUE)
    expect_equal(class(p2$layers[[1]]$geom)[1], 'GeomCol')    
})
