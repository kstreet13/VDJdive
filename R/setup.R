
#' @title Load 10X CellRanger V(D)J data
#' @name readVDJcontigs
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
#' @import S4Vectors
#' @import IRanges
#' @export
setMethod(f = 'readVDJcontigs',
          signature = signature(samples = "character"),
          definition = function(samples, sample.names = names(samples)){
              if(is.null(sample.names)){
                  sample.names <- basename(samples)
              }
              # get TCR contigs
              contigs <- NULL
              for (ii in 1:length(samples)) {
                  tcr_ii <- read.csv(file.path(samples[ii],
                                               "filtered_contig_annotations.csv"))
                  tcr_ii$barcode <- gsub("1$", ii, tcr_ii$barcode)
                  #tcr_ii$type <- "TCR"
                  #tcr_ii$clonotype_tcr <- paste0(sample_ids[ii], "-", tcr_ii$raw_clonotype_id)
                  tcr_ii$sample <- sample.names[ii]
                  contigs <- rbind(contigs, tcr_ii)
              }
              rm(ii, tcr_ii)
              # convert to SplitDataFrameList
              tcr.list <- split(DataFrame(contigs), factor(contigs$barcode))
              
              return(tcr.list)
          })




#' @title Add 10X CellRanger V(D)J data to SingleCellExperiment
#' @name addVDJtoSCE
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
#' @import S4Vectors
#' @importFrom SummarizedExperiment colData
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
              message(loss, ' cells with V(D)J data were dropped because',
                      'they had no match in SingleCellExperiment object (',
                      format(100 * pct, digits = 2),
                      '% of V(D)J data).')
              
              sce$contigs <- tcr.list
              return(sce)
          })
          
#' @rdname addVDJtoSCE
#' @export          
setMethod(f = 'addVDJtoSCE',
          signature = signature(samples = "character", 
                                sce = "SingleCellExperiment"),
          definition = function(samples, sce, sample.names = names(samples), barcode = 'Barcode'){
              # get TCR contigs
              contigs <- readVDJcontigs(samples, sample.names = sample.names)
              return(addVDJtoSCE(contigs, sce, sample.names = sample.names, barcode = barcode))
          })


