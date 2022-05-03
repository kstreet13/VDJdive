#' @title Create a box plot for diversity measures
#' @param ... additional arguments.
#' @name boxVDJ
#' @export
setGeneric(name = "boxVDJ",
           signature = "d",
           def = function(d, ...) standardGeneric("boxVDJ"))

#' @rdname boxVDJ
#' 
#' @description \code{boxVDJ} creates a box plot of the specified diversity. 
#' 
#' @param d A \code{matrix} created with \code{calculateDiversity}.
#' @param sampleGroups A \code{matrix} or \code{data.frame} that 
#' identifies the groups that each sample belongs to. The matrix must contain 
#' two columns. The first column lists the individual samples and should be 
#' called "Sample". The second column should list the group that each sample
#' belongs to (e.g. Normal and Tumor) and be called "Group". If no 
#' sampleGroups dataset is provided, all of the samples will be plotted
#' in one group. 
#' @param method Identifies the type of diversity that is to be plotted. 
#' @param title Character vector with an optional title. 
#' @param legend If TRUE, a legend will be included with the plot. If FALSE,
#' no legend is included in the plot.
#' 
#' @return Returns a \code{ggplot} plot with a box plot that shows the
#' diversity for each sample. A box plot is created for each of the 
#' grouping variables. The individual diversity measures are 
#' plotted on the box plots. 
#'
#' @examples
#' data(contigs)
#' x <- clonoStats(contigs)
#' d <- calculateDiversity(x)
#' sampleGroups <- data.frame(Sample = c("sample1", "sample2"), 
#'                            Group = c("Cancer", "Normal"))
#' boxVDJ(d, sampleGroups = sampleGroups, method = "shannon", 
#'        title = "Shannon diversity", legend = FALSE)
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter scale_color_manual theme_bw labs
#' @export
setMethod(f = "boxVDJ",
          signature = signature(d = "matrix"),
          definition = function(d, sampleGroups = NULL, 
                                method = c("shannon", "simpson", "invsimpson", 
                                           "chao1", "chaobunge"), 
                                title = NULL, legend = FALSE){
              if (!is.null(sampleGroups)) {
                  sampleGroups <- data.frame(sampleGroups)
                  colnames(sampleGroups) <- tolower(colnames(sampleGroups))
              }
              
              if (method == "chaobunge") method <- "chaobunge.estimate"
              method <- paste0("^", method)
              
              if (!is.null(sampleGroups)) {
                  d1 <- t(d[grepl(method, dimnames(d)[[1]], 
                                  ignore.case = TRUE), ])
                  sampleGroups$Diversity <- unlist(d1[match(dimnames(d1)[[2]],
                                                        sampleGroups$sample)])
              } else {
                  d1 <- t(d[grepl(method, dimnames(d)[[1]], 
                                  ignore.case = TRUE), ])
                  sampleGroups <- data.frame(Diversity = unlist(d1), 
                                             group = " ")
              }
              
              cols <- RColorBrewer::brewer.pal(8, name = "Dark2")[seq_along(
                  unique(sampleGroups$group))]
              
              g <- ggplot2::ggplot(sampleGroups, ggplot2::aes(x = group, 
                                                              y = Diversity, 
                                                              color = group)) 
              
              if (legend) {
                  g <- g +
                      ggplot2::geom_boxplot() +
                      ggplot2::geom_jitter(height = 0, width = 0.1) +
                      ggplot2::scale_color_manual(name = "", values = cols) +
                      ggplot2::theme_bw()
                  
              } else {
                  g <- g + 
                      ggplot2::geom_boxplot(show.legend = FALSE) +
                      ggplot2::geom_jitter(height = 0, width = 0.1, 
                                           show.legend = FALSE) +
                      ggplot2::scale_color_manual(name = "", values = cols) +
                      ggplot2::theme_bw()
              }
              
              if (!is.null(title)) {
                  g <- g + ggplot2::labs(title = title, x = NULL)
              } else {
                  g <- g + ggplot2::labs(x = NULL)
              }
              g
          })

