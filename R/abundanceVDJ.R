#' @include clonoStats_class.R
NULL

#' @title Create an abundance graph for clonotype expansion
#' @param ... additional arguments.
#' @name abundanceVDJ
#' @export
setGeneric(name = "abundanceVDJ",
           signature = "x",
           def = function(x, ...) standardGeneric("abundanceVDJ"))


#' @rdname abundanceVDJ
#'
#' @description \code{abundanceVDJ} creates a dot plot using ggplot that shows the
#'   number of reads for each clonotype in each sample and labels the most 
#'   abundant clonotypes.
#' 
#' @param x A \code{matrix} created with \code{clonoStats}.
#' @param annotate An integer that specifies how many of the most abundant
#' clonotypes should be annotated on the plot. 
#' @param title Character vector with an optional title. If FALSE, no title
#' is generated.
#' 
#' @return Returns a \code{ggplot} plot with a dot plot that shows the
#' abundance of the clonotypes in each sample. The most abundant 
#' clonotypes are annotated on the plot and ordered from most abundant 
#' to least abundant.
#'
#' @examples
#' data('contigs')
#' x <- clonoStats(contigs)
#' abundanceVDJ(x)
#' 
#' @importFrom ggplot2 ggplot geom_text aes labs theme_bw geom_jitter element_blank
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom cowplot panel_border
#' @importFrom gridExtra grid.arrange
#' @export
setMethod(f = "abundanceVDJ",
          signature = signature(x = "clonoStats"),
          definition = function(x, annotate = 5, title = NULL) {
              a <- clonoAbundance(x)
              rownames(a) <- clonoNames(x)
              dat <- NULL
              dat2 <- NULL
              nms <- colnames(a)
              for (i in seq_len(ncol(a))) {
                  tmp <- a[a[, i] > 0.05, nms[i]]
                  tmp <- c(sort(tmp, decreasing = TRUE))
                  dat <- rbind(dat, data.frame(clonotype = seq_along(tmp), count = tmp, 
                                               sample = nms[i]))
                  len <- min(annotate, length(tmp))
                  dat2 <- rbind(dat2, data.frame(clono = c(" ", names(tmp)[seq(1, len)], ""), 
                                                 sample = nms[i], 
                                                 y = seq(len + 2, 1), 
                                                 ylab = c(" ", seq(1, len), " ")))
              }
              
              g <- ggplot(dat, aes(x = sample, y = count)) +
                  geom_jitter(height = 0, width = 0.1, color = "#1b9e77", alpha = 0.5) +
                  labs(title = NULL, x = NULL, y = "Abundance") +
                  theme_bw()
              
              g2 <- ggplot(dat2, aes(x = sample, y = y)) +
                  geom_text(label = dat2$clono) +
                  labs(title = title, x = NULL, y = NULL) +
                  theme_bw() +
                  panel_border(remove = TRUE) +
                  theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), 
                                 axis.line = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank()
                  ) 
                  
              grid.arrange(g2, g, layout_matrix = matrix(c(1, 2, 2, 2), ncol = 1))
              
          })

