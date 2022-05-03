#' @include clonoStats_class.R
NULL

#' @title Create a bar graph for clonotype expansion
#' @param ... additional arguments.
#' @name barVDJ
#' @export
setGeneric(name = "barVDJ",
           signature = "x",
           def = function(x, ...) standardGeneric("barVDJ"))


#' @rdname barVDJ
#'
#' @description \code{barVDJ} creates a barplot using ggplot that shows the
#'   number of reads in the sample and colors the sample in accordance to the
#'   amount of diversity.
#' 
#' @param x A \code{matrix} created with \code{summarizeClonotypes}.
#' @param bySample If TRUE, the plot will be separated by sample. If FALSE,
#' all samples will be combined as if they were one sample.
#' @param title Character vector with an optional title. If FALSE, no title
#' is generated.
#' @param legend If TRUE, a legend will be included with the plot. If FALSE,
#' no legend is included in the plot.
#' 
#' @return Returns a \code{ggplot} plot with a barplot that shows the
#' abundance of the clonotypes. The coloring indicates the number of cells
#' for each clonotype with darker colors being clonotypes with a single cell
#' (singletons) and lighter colors having more cells with that clonotype
#' (expanded clonotype).
#'
#' @examples
#' data('contigs')
#' x <- clonoStats(contigs)
#' barVDJ(x)
#' 
#' @importFrom ggplot2 ggplot geom_col aes scale_fill_continuous labs theme_bw
#' @importClassesFrom Matrix dgCMatrix
#' @export
setMethod(f = "barVDJ",
          signature = signature(x = "Matrix"),
          definition = function(x, bySample = TRUE,
                                title = NULL, legend = FALSE){
              dat <- NULL
              nms <- colnames(x)
              for (i in seq_len(ncol(x))) {
                  tmp <- x[x[, i] > 0.5, nms[i]]
                  tmp <- c(sort(tmp, decreasing = TRUE),0)
                  dat <- rbind(dat, data.frame(clonotype = seq_along(tmp), count = tmp, 
                                               sample = nms[i]))
              }
              g <- ggplot2::ggplot(dat, ggplot2::aes(x = sample, y = count, 
                                                     weight = count, fill = log(count)))
              
              if (legend) {
                  g <- g + ggplot2::geom_col(position = "stack")
              } else {
                  g <- g + ggplot2::geom_col(position = "stack", show.legend = FALSE)
              }
              
              g <- g +
                  ggplot2::scale_fill_continuous(type = "viridis") +
                  ggplot2::labs(x = NULL, y = "Number of T cells", title = title) +
                  ggplot2::theme_bw()
              g
          })

#' @rdname barVDJ
#' @export
setMethod(f = "barVDJ",
          signature = signature(x = "matrix"),
          definition = function(x, ...){
              barVDJ(Matrix(x), ...)
          })

#' @rdname barVDJ
#' @export
setMethod(f = "barVDJ",
          signature = signature(x = "clonoStats"),
          definition = function(x, ...){
              barVDJ(x@abundance, ...)
          })
