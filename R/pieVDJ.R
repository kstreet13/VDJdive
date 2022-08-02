#' @include clonoStats_class.R
NULL

#' @title Create a pie chart for clonotype expansion
#' @param ... additional arguments.
#' @name pieVDJ
#' @export
setGeneric(name = "pieVDJ",
           signature = "x",
           def = function(x, ...) standardGeneric("pieVDJ"))



#' @rdname pieVDJ
#'
#' @description \code{pieVDJ} creates a list of pie charts created using ggplot
#'   that shows the the level of expansion in each clonotype.
#' 
#' @param x A \code{matrix} created with \code{summarizeClonotypes}.
#' @param legend Can take on the values use in the legend.position command
#' in ggplot ("left","top", "right", "bottom", or a numeric vector) to 
#' indicate where the legend should be placed. If left NULL, no legend will 
#' be created.
#' 
#' @return If \code{x} contains more than one sample, a list of pie charts
#' will be returned. If \code{x} contains only one sample, the pie chart
#' will be returned. The coloring indicates the number of cells
#' for each clonotype with darker colors being clonotypes with a single cell
#' (singletons) and lighter colors having more cells with that clonotype
#' (expanded clonotype).
#'
#' @examples
#' data('contigs')
#' x <- clonoStats(contigs)
#' pieVDJ(x)
#' 
#' @importFrom ggplot2 ggplot aes geom_col coord_polar scale_fill_continuous labs theme_void theme
#' @importClassesFrom Matrix dgCMatrix
#' @export
setMethod(f = "pieVDJ",
          signature = signature(x = "Matrix"),
          definition = function(x, legend = "bottom"){
              pieG <- NULL
              
              for (i in seq_len(ncol(x))) {
                  numcells <- sum(x[, i], na.rm = TRUE)
                  singletons <- c(sum(x[x[, i] <= 1, i]), 1)
                  nm <- colnames(x)[i]
                  if (!(numcells - singletons[1])) {
                      dat <- data.frame(count = 1, col = 1, clonotype = 1)
                      warning("Sample ", nm, " only contains singletons")
                  } else {
                      cl <- x[x[, i] > 1, i]
                      dat <- data.frame(count = cl[order(cl)])
                      dat$col <- dat$count
                      dat <- rbind(singletons, dat)
                      dat$clonotype <- seq_len(nrow(dat))
                  }
                  
                  pieG[[i]] <- ggplot(dat, aes(x = "", y = count, fill = log(col), weight = count))
                  
                  if (!is.null(legend)) {
                      pieG[[i]] <- pieG[[i]] + geom_col(position = "stack")
                  } else {
                      pieG[[i]] <- pieG[[i]] + geom_col(position = "stack", show.legend = FALSE)
                  }
                  
                  pieG[[i]] <- pieG[[i]] + 
                      coord_polar("y", start = pi / 2) + 
                      scale_fill_continuous("", 
                                                     type = "viridis", 
                                                     breaks = c(0, max(log(dat$col))), 
                                                     labels = c("Singletons", "Most expanded")) +
                      labs(x = NULL, y = NULL, fill = NULL, title = nm) +
                      theme_void() +
                      theme(legend.position = legend)
                  
              }
              
              if (length(pieG) == 1) {
                  pieG[[1]]
              } else {
                  pieG
              }
          })


#' @rdname pieVDJ
#' @export
setMethod(f = "pieVDJ",
          signature = signature(x = "matrix"),
          definition = function(x, ...){
              pieVDJ(Matrix(x), ...)
          })

#' @rdname pieVDJ
#' @export
setMethod(f = "pieVDJ",
          signature = signature(x = "clonoStats"),
          definition = function(x, ...){
              pieVDJ(clonoAbundance(x), ...)
          })

