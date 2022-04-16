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
#' @return Returns a \code{ggplot} plot with a box plot that shows the
#' diversity for each sample. A box plot is created for each of the 
#' grouping variables. The individual diversity measures are 
#' plotted on the box plots. 
#'
#' @examples
#' data(contigs)
#' samples <- vapply(contigs[,'sample'], function(x){ x[1] }, 'A')
#' counts <- EMquant(contigs)
#' x <- t(summarizeClonotypes(counts, samples))
#' d <- calculateDiversity(x)
#' 
#' sampleGroups <- data.frame(Sample = c("sample1", "sample2"), 
#'                            Group = c("Cancer", "Normal"))
#' boxVDJ(d, sampleGroups = sampleGroups, method = "shannon", 
#'        title = "Shannon diversity", legend = FALSE)
#' @export
#'

divLineVDJ <- function(d, 
                       sampleGroups = NULL, 
                       method = c("shannon", "simpson", "invsimpson", "chao1", "chaobunge"), 
                       title = NULL) {
  
  if (method == "chaobonge") {
    # method <- "chaobonge.estimate"
    # error <- "chaobonge.error"
    dat <- d[grepl(method, rownames(d), ignore.case = TRUE), ]
    dat <- data.frame(type = gsub("chaobonge.", "", rownames(dat)), dat)
    dat <- tidyr::pivot_longer(dat, cols = colnames(dat)[-1], 
                               names_to = "Sample", 
                               values_to = "Diversity") %>%
      tidyr::pivot_wider(names_from = type, 
                         values_from = Diversity)
    colnames(dat)[colnames(dat) == "estimate"] <- "Diversity"
    dat$Diversity <- dat$Diversity / 100
    dat$error <- dat$error / 100
    dat$pairs <- gsub(normal, "", dat$Sample)
    dat$pairs <- gsub(tumor, "", dat$pairs)
    dat$type <- ifelse(grepl(tumor, dat$Sample), "Tumor", "Normal")
    dat$position <- ifelse(dat$type == "Tumor", 2, 1)
    samples <- unique(dat$pairs)
    count <- 0
    for (i in 1:length(samples)) {
      dat$position[dat$pairs == samples[i]] <- dat$position + count
      count <- count + 0.1
    }
    
    ggplot(dat, aes(x = position, y = Diversity, color = type, group = pairs)) +
      geom_line(color = "gray40") +
      geom_point() +
      geom_errorbar(aes(ymin = Diversity - error, ymax = Diversity + error), 
                    width = 0.04, position = position_dodge2(padding = 3),
                    size = 1) +
      scale_color_manual(name = "", values = c("blue", "red")) +
      labs(x = NULL, title = title) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) 
    
  } else {
    m <- paste0("^", method)
    dat <- unlist(d[grepl(m, rownames(d), ignore.case = TRUE), ])
    dat <- tidyr::pivot_longer(dat, everything(), 
                               names_to = "Sample", 
                               values_to = "Diversity") 
    dat$pairs <- gsub(normal, "", dat$Sample)
    dat$pairs <- gsub(tumor, "", dat$pairs)
    dat$type <- ifelse(grepl(tumor, dat$Sample), "Tumor", "Normal")
    
    ggplot(dat, aes(x = type, y = Diversity, color = type)) +
      geom_line(aes(group = pairs), color = "gray40") +
      geom_point() +
      scale_color_manual(name = "", values = c("blue", "red")) +
      labs(x = NULL, title = title) +
      theme_bw()
  }
}
