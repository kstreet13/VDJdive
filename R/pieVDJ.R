#' pieVDJ creates a pie chart for clonotype expansion
#'
#' \code{pieVDJ} creates a list of pie charts created using ggplot that 
#' shows the the level of epxansion in each clonotype.
#' @param x A \code{matrix} created with \code{summarizeClonotypes}.
#' @param legend Can take on the values use in the legend.position command
#' in ggplot ("left","top", "right", "bottom", or a numeric vector) to 
#' indicate where the legend should be placed. If left NULL, no legend will 
#' be created. .
#' @return If \code{x} contains more than one sample, a list of pie charts
#' will be returned. If \code{x} contains only one sample, the pie chart
#' will be returned. The coloring indicates the number of cells
#' for each clonotype with darker colors being clonotypes with a single cell
#' (singletons) and lighter colors having more cells with that clonotype
#' (expanded clonotype).
#'
#' @examples
#' data(contigs)
#' samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
#' counts <- EMquant(contigs)
#' x <- t(summarizeClonotypes(counts, samples))
#' pieVDJ(x)
#' @export
#'
pieVDJ <- function(x, legend = "bottom") {
  
  pieG <- NULL
  
  for (i in 1:ncol(x)) {
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
      dat$clonotype <- 1:nrow(dat)
    }
    
    pieG[[i]] <- ggplot2::ggplot(dat, ggplot2::aes(x = "", y = count, fill = log(col), weight = count))
    
    if (!is.null(legend)) {
      pieG[[i]] <- pieG[[i]] + ggplot2::geom_col(position = "stack")
    } else {
      pieG[[i]] <- pieG[[i]] + ggplot2::geom_col(position = "stack", show.legend = FALSE)
    }
    
    pieG[[i]] <- pieG[[i]] + 
      ggplot2::coord_polar("y", start = pi / 2) + 
      ggplot2::scale_fill_continuous("", 
                                     type = "viridis", 
                                     breaks = c(0, max(log(dat$col))), 
                                     labels = c("Singletons", "Most expanded")) +
      ggplot2::labs(x = NULL, y = NULL, fill = NULL, title = nm) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = legend)
    
  }
  
  # Make legend
  # datL <- data.frame(values = 1:100)
  # pieL <- ggplot2::ggplot(datL, ggplot2::aes(x = "", y = values, fill = values, weight = values)) +
  #   ggplot2::geom_col(position = "stack") +
  #   ggplot2::scale_fill_continuous("", type = "viridis", breaks = c(1, 100), 
  #                                  labels = c("Singletons", "Most expanded")) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(legend.position = "bottom") 
  

  if (length(pieG) == 1) {
    pieG[[1]]
  } else {
    pieG
  }
}

