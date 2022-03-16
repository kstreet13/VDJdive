
#' @title Load 10X CellRanger V(D)J data
#' @name readVDJcontigs
#' @param ... additional arguments.
#' @export
setGeneric(name = "readVDJcontigs",
           signature = "samples",
           def = function(samples, ...) standardGeneric("readVDJcontigs"))

#' @rdname readVDJcontigs
#'
#' @description Creates a \code{SplitDataFrameList} (see
#'   \code{\link[IRanges]{DataFrameList}}) from a vector of directory names
#'   corresponding to the output of the CellRanger V(D)J pipeline.
#'
#' @param samples A character vector containing one or more directory names,
#'   each corresponding to a 10X sample. Each directory should contain a file
#'   named \code{filtered_contig_annotations.csv}.
#' @param sample.names A character vector of length equal to \code{samples}
#'   containing the sample names to store in the output object. If \code{NULL},
#'   the \code{basename}s of each directory will be used.
#'
#' @details The resulting list of DataFrames contains all the data in
#'   \code{filtered_contig_annotations.csv}, split by cell barcode. Note that
#'   the index of each sample in \code{samples} is concatenated to the cell
#'   barcodes, so that cells from different samples cannot have identical
#'   barcodes.
#'
#' @return A \code{SplitDataFrameList} object containing data on all contigs,
#'   grouped by cell barcode.
#'
#' @examples
#' # write the example data to a temporary directory
#' example(writeVDJcontigs)
#' # specify sample locations and read in data
#' samples <- file.path(loc, c('sample1','sample2'))
#' contigs <- readVDJcontigs(samples)
#'
#' @importFrom utils read.csv
#' @import S4Vectors
#' @importClassesFrom IRanges SplitDataFrameList
#' @export
setMethod(f = 'readVDJcontigs',
          signature = signature(samples = "character"),
          definition = function(samples, sample.names = names(samples)){
              if(is.null(sample.names)){
                  sample.names <- basename(samples)
              }
              # get TCR contigs
              contigs <- NULL
              for (ii in seq_along(samples)) {
                  tcr_ii <- read.csv(file.path(
                      samples[ii], "filtered_contig_annotations.csv"))
                  tcr_ii$barcode <- gsub("1$", ii, tcr_ii$barcode)
                  #tcr_ii$type <- "TCR"
                  #tcr_ii$clonotype_tcr <- paste0(sample_ids[ii], "-", 
                  #tcr_ii$raw_clonotype_id)
                  tcr_ii$sample <- sample.names[ii]
                  contigs <- rbind(contigs, tcr_ii)
              }
              rm(ii, tcr_ii)
              # convert to logical ("None" treated as FALSE)
              contigs$is_cell <- contigs$is_cell %in% c('true','True','TRUE')
              contigs$high_confidence <- contigs$high_confidence %in% 
                  c('true','True','TRUE')
              contigs$full_length <- contigs$full_length %in% 
                  c('true','True','TRUE')
              contigs$productive <- contigs$productive %in% 
                  c('true','True','TRUE')
              # convert to SplitDataFrameList
              tcr.list <- split(DataFrame(contigs), factor(contigs$barcode))

              return(tcr.list)
          })




#' @title Add 10X CellRanger V(D)J data to SingleCellExperiment
#' @name addVDJtoSCE
#' @param ... additional arguments.
#' @export
setGeneric(name = "addVDJtoSCE",
           signature = c("samples","sce"),
           def = function(samples, sce, ...) standardGeneric("addVDJtoSCE"))

#' @rdname addVDJtoSCE
#'
#' @description Matches CellRanger V(D)J data to paired data in an existing
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#'
#' @param samples A character vector containing one or more directory names,
#'   each corresponding to a 10X sample. Each directory should contain a file
#'   named \code{filtered_contig_annotations.csv}. Alternatively, a
#'   \code{SplitDataFrameList}, the output of \code{\link{readVDJcontigs}}.
#' @param sce A \code{SingleCellExperiment} object.
#' @param sample.names A character vector of length equal to \code{samples}
#'   containing the sample names to store in the output object. If \code{NULL}
#'   and \code{samples} is a character vector, the \code{basename}s of each
#'   directory will be used.
#' @param barcode The column name from the \code{colData} of \code{sce}
#'   containing cell barcodes. These should match the barcodes in the V(D)J data
#'   (see Details). Alternatively, a vector of cell barcodes of length equal to
#'   \code{ncol(sce)}.
#'
#' @details Matching cell barcodes between data objects can cause problems,
#'   because different methods have different ways of ensuring barcodes are
#'   unique across all samples. This function and \code{\link{readVDJcontigs}}
#'   follow the naming conventions of \code{\link[DropletUtils]{read10xCounts}},
#'   where the sample index (in the \code{samples} vector) is appended to each
#'   cell barcode, to ensure that each barcode is unique, across all samples. If
#'   \code{sce} was created by a different method, such as conversion from a
#'   \code{Seurat} object, you may need to check the barcode naming convention.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'   with an element named \code{contigs} added to the \code{colData},
#'   representing the V(D)J data.
#'
#' @examples
#' # load example V(D)J data
#' data('contigs')
#' # make SCE object with matching barcodes and sample IDs
#' ncells <- 24
#' u <- matrix(rpois(1000 * ncells, 5), ncol = ncells)
#' barcodes <- vapply(contigs[,'barcode'], function(x){ x[1] }, 'A')
#' samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = u),
#'                             colData = data.frame(Barcode = barcodes,
#'                                                  sample = samples))
#' sce <- addVDJtoSCE(contigs, sce)
#' sce$contigs
#'
#' @import S4Vectors
#' @importFrom SummarizedExperiment colData colData<-
#' @importClassesFrom IRanges SplitDataFrameList
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
setMethod(f = 'addVDJtoSCE',
          signature = signature(samples = "SplitDataFrameList",
                                sce = "SingleCellExperiment"),
          definition = function(samples, sce, sample.names = NULL,
                                barcode = 'Barcode'){
              contigs <- unlist(samples)
              if(length(barcode) == 1){
                  stopifnot(barcode %in% names(colData(sce)))
                  bcVar <- as.character(sce[[barcode]])
              }else{
                  stopifnot(length(barcode) == ncol(sce))
                  bcVar <- as.character(barcode)
              }

              tcr.list <- split(DataFrame(contigs),
                                factor(contigs$barcode, bcVar))
              loss <- length(unique(contigs$barcode)) -
                  sum(lengths(tcr.list)>0)
              pct <- loss / length(unique(contigs$barcode))
              if(loss > 0){
                  message(loss, ' cells with V(D)J data were dropped because ',
                          'they had no match in SingleCellExperiment object (',
                          format(100 * pct, digits = 2),
                          '% of V(D)J data).')
              }
              sce$contigs <- tcr.list
              return(sce)
          })

#' @rdname addVDJtoSCE
#' @export
setMethod(f = 'addVDJtoSCE',
          signature = signature(samples = "character",
                                sce = "SingleCellExperiment"),
          definition = function(samples, sce, sample.names = names(samples), 
                                barcode = 'Barcode'){
              # get TCR contigs
              contigs <- readVDJcontigs(samples, sample.names = sample.names)
              return(addVDJtoSCE(contigs, sce, sample.names = sample.names, 
                                 barcode = barcode))
          })








#' @title Write V(D)J contig data in 10X format
#' @name writeVDJcontigs
#' @param ... additional arguments.
#' @export
setGeneric(name = "writeVDJcontigs",
           signature = c("path","x"),
           def = function(path, x, ...) standardGeneric("writeVDJcontigs"))


#' @rdname writeVDJcontigs
#'
#' @description Write V(D)J data to a series of directories, each containing a
#'   CSV file with the data for an individual sample.
#'
#' @param path A string containing the path to the output directory.
#' @param x A \code{SplitDataFrameList} object containing V(D)J contig
#'   information, split by cell barcodes, as created by \code{readVDJcontigs}.
#'
#' @returns
#' Creates various subdirectories of the directory specified in \code{path}.
#' Each subdirectory is named for a sample found in the dataset, \code{x}, and
#' contains a CSV filed named \code{filtered_contig_annotations.csv}.
#'
#' @examples
#' data('contigs')
#' loc <- tempdir()
#' writeVDJcontigs(loc, contigs)
#'
#' @importFrom utils write.csv
#' @export
setMethod(f = 'writeVDJcontigs',
          signature = signature(path = "character", x = "SplitDataFrameList"),
          definition = function(path, x){

    contigs <- data.frame(unlist(x))

    if(!dir.exists(path)){
        dir.create(path)
    }
    for(samp in unique(contigs$sample)){
        contigs.ii <- contigs[which(contigs$sample==samp), ]
        contigs.ii$sample <- NULL
        path.ii <- file.path(path, samp)
        if(!dir.exists(path.ii)){
            dir.create(path.ii)
        }
        write.csv(contigs.ii, file.path(path.ii,
                                        'filtered_contig_annotations.csv'))
    }
})



