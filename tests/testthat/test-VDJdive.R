context("Test VDJdive.")
library(utils) # needed for data()
library(stats) # for rpois()

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
                             sample = samples))
    sce <- addVDJtoSCE(contigs, sce)
    sce <- uniqueQuant(sce)
    
    # splitClonotypes
    countsList <- splitClonotypes(sce, by = 'sample')
    expect_equal(length(countsList), 2)
    expect_equivalent(dim(countsList[[1]]), c(12, 7))
    expect_equivalent(dim(countsList[[2]]), c(12, 7))
    expect_equal(sum(countsList[[1]]), 6)
    
    # summarizeClonotypes
    sampleLevelCounts <- summarizeClonotypes(sce, by = 'sample')
    expect_equivalent(dim(sampleLevelCounts), c(7, 2))
    expect_equivalent(colSums(sampleLevelCounts), c(6, 6))
    
    sampleLevelCounts <- summarizeClonotypes(sce, by = 'sample', mode = 'tab')
    expect_equivalent(dim(sampleLevelCounts), c(5, 2))
    expect_equivalent(colSums(sampleLevelCounts), c(7, 7))
    expect_equivalent(rownames(sampleLevelCounts), as.character(0:4))
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
                             sample = samples))

    sce <- addVDJtoSCE(contigs, sce)
    expect_true('contigs' %in% names(colData(sce)))
    expect_equivalent(lengths(sce$contigs),
                      c(3,6,2,4,2,3,2,2,2,4,2,2,3,1,2,4,2,2,1,2,2,1,1,3))
})

test_that("assignment functions work", {
    # load example data
    data("contigs")

    Uniqcounts <- uniqueQuant(contigs)
    expect_equivalent(dim(Uniqcounts), c(24, 7))
    expect_equal(sum(Uniqcounts), 12)

    CRcounts <- CRquant(contigs)
    expect_equivalent(dim(CRcounts), c(24, 19))
    expect_equal(sum(CRcounts), 22)

    EMcounts <- EMquant(contigs)
    expect_equivalent(dim(EMcounts), c(24, 60))
    expect_equal(sum(EMcounts), 22)
    
    EMcounts2 <- EMquant(contigs, method = 'r')
    expect_equivalent(dim(EMcounts2), c(24, 60))
    expect_equal(sum(EMcounts2), 22)
    expect_true(max(abs(EMcounts2 - EMcounts)) < .0001)

    # make SCE object with matching barcodes and sample IDs
    ncells <- 24
    u <- matrix(rpois(1000 * ncells, 5), ncol = ncells)
    barcodes <- vapply(contigs[,'barcode'], function(x){ x[1] }, 'A')
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = u),
        colData = data.frame(Barcode = barcodes,
                             sample = samples))
    sce <- addVDJtoSCE(contigs, sce)
    
    sceUniq <- uniqueQuant(sce)
    expect_true(max(abs(sceUniq$clono - Uniqcounts)) < .0001)
    
    sceCR <- CRquant(sce)
    expect_true(max(abs(sceCR$clono - CRcounts)) < .0001)
    
    sceEM <- EMquant(sce)
    expect_true(max(abs(sceEM$clono - EMcounts)) < .0001)
    
    sceEMsamp <- EMquant(sce, sample = 'sample')
    sceEMsamp2 <- EMquant(sce, sample = sce$sample)
    expect_true(max(abs(sceEMsamp$clono - sceEMsamp2$clono)) < .0001)
})

test_that("assignment functions handle edge cases", {
    # load example data
    data("contigs")
    
    # empty sample
    sample <- factor(vapply(contigs[,'sample'], function(x){ x[1] }, 'A'))
    levels(sample) <- c('sample1','sample2','sample3')
    Uniqcounts <- uniqueQuant(contigs, sample = sample)
    expect_equivalent(dim(Uniqcounts), c(24, 7))
    expect_equal(sum(Uniqcounts), 12)
    
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
                             sample = samples))
    sce <- addVDJtoSCE(contigs, sce)
    sceUniq <- uniqueQuant(sce, sample ='sample')
    expect_equivalent(dim(sceUniq$clono), c(30, 7))
    expect_equal(sum(sceUniq$clono), 12)
    expect_equal(sum(sceUniq$clono[sce$sample=='sample3',]), 0)
    
    # BCR data and mixtures
    Uniqcounts <- uniqueQuant(contigs, type = 'BCR')
    expect_equivalent(dim(Uniqcounts), c(24, 0))
    contigs[contigs[,'chain']=='TRA','chain'] <- 'IGH'
    Uniqcounts <- uniqueQuant(contigs, type = 'BCR')
    expect_equivalent(dim(Uniqcounts), c(24, 0))
    contigs[contigs[,'chain']=='TRB','chain'] <- 'IGL'
    
    # let it detect BCR
    Uniqcounts <- uniqueQuant(contigs)
    expect_equivalent(dim(Uniqcounts), c(24, 7))
    expect_equal(sum(Uniqcounts), 12)
    
    # reset
    rm(contigs)
    data("contigs")
    
    # empty sample
    EMcounts <- EMquant(contigs, sample = sample)
    expect_equivalent(dim(EMcounts), c(24, 60))
    expect_equal(sum(EMcounts), 22)
    
    # SCE object with EXTRA CELLS and SAMPLE
    sceEM <- EMquant(sce, sample ='sample')
    expect_equivalent(dim(sceEM$clono), c(30, 60))
    expect_equal(sum(sceEM$clono), 22)
    expect_equal(sum(sceEM$clono[sce$sample=='sample3',]), 0)
    
    # BCR data and mixtures
    EMcounts <- EMquant(contigs, type = 'BCR')
    expect_equivalent(dim(EMcounts), c(24, 0))
    contigs[contigs[,'chain']=='TRA','chain'] <- 'IGH'
    EMcounts <- EMquant(contigs, type = 'BCR')
    expect_equivalent(dim(EMcounts), c(24, 0))
    contigs[contigs[,'chain']=='TRB','chain'] <- 'IGL'
    
    # let it detect BCR
    EMcounts <- EMquant(contigs)
    expect_equivalent(dim(EMcounts), c(24, 60))
    expect_equal(sum(EMcounts), 22)
})

test_that("diversity calculation works", {
    # load example data
    data("contigs")
    CRcounts <- CRquant(contigs)
    EMcounts <- EMquant(contigs)
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    
    kCR <- summarizeClonotypes(CRcounts, by = samples)
    kEM <- summarizeClonotypes(EMcounts, by = samples)
    
    divCR <- calculateDiversity(kCR)
    expect_equivalent(dim(divCR), c(6, 2))
    expect_equivalent(colnames(divCR), c('sample1', 'sample2'))

    expect_warning({
        divEM <- calculateDiversity(kEM)
    }, regexp = 'Cut-off was too low')
    expect_equivalent(dim(divEM), c(6, 2))
    expect_equivalent(colnames(divEM), c('sample1', 'sample2'))
})

test_that("plotting functions work", {
    # load example data
    data("contigs")
    
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    counts <- EMquant(contigs)
    x <- summarizeClonotypes(counts, samples)
    p1 <- barVDJ(x)
    expect_equal(class(p1$layers[[1]]$geom)[1], 'GeomCol')
    
    p2 <- barVDJ(x, bySample = FALSE, title = 'bar plot', legend = TRUE)
    expect_equal(class(p2$layers[[1]]$geom)[1], 'GeomCol')    
})
