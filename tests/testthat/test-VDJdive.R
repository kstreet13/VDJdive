context("Test VDJdive.")
library(utils) # needed for data()
library(stats) # for rpois() and runif()

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

test_that("quantification functions work", {
    # load example data
    data("contigs")

    counts <- uniqueQuant(contigs)
    expect_equivalent(dim(counts), c(24, 7))
    expect_equal(sum(counts), 12)

    counts <- CRquant(contigs)
    expect_equivalent(dim(counts), c(24, 19))
    expect_equal(sum(counts), 22)

    counts <- EMquant(contigs)
    expect_equivalent(dim(counts), c(24, 60))
    expect_equal(sum(counts), 22)

    counts2 <- EMquant(contigs, method = 'r')
    expect_equivalent(dim(counts2), c(24, 60))
    expect_equal(sum(counts2), 22)
    expect_true(max(abs(counts2 - counts)) < .001)

})

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
    expect_equivalent(dim(sampleLevelCounts), c(2, 7))
    expect_equivalent(rowSums(sampleLevelCounts), c(6, 6))
})

test_that("plotting functions work", {
    # load example data
    data("contigs")
    
    samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
    counts <- EMquant(contigs)
    x <- t(summarizeClonotypes(counts, samples))
    p1 <- barVDJ(x)
    expect_equal(class(p1$layers[[1]]$geom)[1], 'GeomCol')
    
    p2 <- barVDJ(x, bySample = FALSE, title = 'bar plot', legend = TRUE)
    expect_equal(class(p2$layers[[1]]$geom)[1], 'GeomCol')    
})
