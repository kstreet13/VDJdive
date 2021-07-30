library(dplyr)
library(vegan)
library(reshape2)
library(fossil)
library(breakaway)

#Start from the tcr_file_fin file

create_counts_matrix <- function(x) {
  # collapse counts by clonotype (i.e., sum by sample and cdr3)
  collapsed <- x %>%
    select(sample, cdr3, reads) %>%
    group_by(sample, cdr3) %>%
    summarise(reads=sum(reads), .groups="drop") %>% # .group silence summarise warning
    as.data.frame()
  
  # melt and make data.frame with one column for each sample
  melted <- melt(collapsed, id.vars=c("sample", "cdr3")) %>%
    dcast(cdr3 ~ sample, value.var = "value", fill=0) #rows cdr3 cols sample value reads
  # set rownames to cdr3 and remove cdr3 column
  rownames(melted) <- melted$cdr3
  melted <- select(melted, -cdr3)
  
  # returns data.frame with cdr3 as rownames, and one column with counts
  # for each sample (colnames are sample names)
  return(melted)
}

#Create matrix k for diversity est
k <- create_counts_matrix(tcr_file_fin)

calculate_diversity <- function(k, methods, ...) { #any parm but k or methods is ... for e.g. cutoff for breakaway which will pass down
  # function to calculate diversity for a given method 
  div_function <- function(k, method, ...) {
    if (method %in% c("simpson", "shannon", "invsimpson")) {
      # simpson, shannon and invsimpson calculated using diversity()
      ret <- data.frame(t(diversity(t(k), ..., index=method)))
      rownames(ret) <- paste0(method, ".estimate") #paste .est after running fun for rownames
    } else if ("chao1" == method) {
      # use fossil::chao1() for chao1
      ret <- data.frame(t(apply(k, 2, fossil::chao1, ...))) #apply chao1 fun for every column
      rownames(ret) <- paste0(method, ".estimate")
    } else if ("chaobonge" == method) {
      # use breakaway() for chaobonge
      ret <- data.frame(apply(k, 2, function(x, ...) {
        bw <- breakaway(x, ...)
        c(estimate=bw$estimate, error=bw$error)
      }))
      rownames(ret) <- paste0(method, ".", rownames(ret))
    } else {
      stop ("invalid diversity method")
    } #error message we can change later
    return(ret)
  }
  
  # loop over methods
  results <- NULL
  for (method in methods) {
    if (is.null(results)) {
      results <- div_function(k, method, ...)#return results of the ret in a dataframe format if it is null
    } else {
      results <- rbind.data.frame(results, div_function(k, method, ...)) #if more than one method rbind all (not null)
    }
  }
  return(results)
}

