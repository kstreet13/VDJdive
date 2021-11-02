#EMDiversity

library(breakaway)
library(fossil)
library(vegan)

FileToy<-readRDS("~/Desktop/PostDoColabProj/toyTCRdata.rds")
quant_Toy<-EMquant(FileToy)



createCountsMatrixEm<-function(x) {
  df<-as.data.frame(colSums(x$clono))
  return(df)
}

k2<-createCountsMatrixEm(quant_Toy)


calculateDiversity <- function(k, methods, ...) {
  # any parm but k or methods is ... for e.g. cutoff for breakaway which will pass down
  # function to calculate diversity for a given method
  
  # Set EM to TRUE of not all values are integers
  EM = !all(as.integer(k[, 1]) == k[, 1])

  div_function <- function(k, method,EM, ...) {
    if (method %in% c("simpson", "shannon", "invsimpson")) {
      # simpson, shannon and invsimpson calculated using diversity() from vegan package
      ret <- data.frame(t(diversity(t(k), ..., index=method)))
      rownames(ret) <- paste0(method, ".estimate") #paste .est after running fun for rownames
    } else if ("chao1" == method & !EM) {
      # use fossil::chao1() for chao1
      ret <- data.frame(t(apply(k, 2, fossil::chao1, ...))) #apply chao1 function for every column
      rownames(ret) <- paste0(method, ".estimate")
    } else if ("chao1" == method  & EM) {
      # use fossil::chao1() for chao1
      ret <- round(k*100,0)
      ret <- data.frame(t(apply(ret, 2, fossil::chao1, ...))) #apply chao1 function for every column
      rownames(ret) <- paste0(method, ".estimate")
    } 
    else if ("chaobonge" == method  & !EM) {
      # use breakaway() for chaobonge
      ret <- data.frame(apply(k, 2, function(x, ...) {
        bw <- breakaway(x, ...)
        c(estimate=bw$estimate, error=bw$error)
      }))
      rownames(ret) <- paste0(method, ".", rownames(ret))
    } 
    else if ("chaobonge" == method & EM) {
      # use breakaway() for chaobonge
      ret <- round(k*100,0)
      ret <- data.frame(apply(ret, 2, function(x, ...) {
        bw <- breakaway(x, ...)
        c(estimate=bw$estimate, error=bw$error)
      }))
      rownames(ret) <- paste0(method, ".", rownames(ret))
    }
    else {
      stop ("invalid diversity method")
    } #error message we can change later
    return(ret)
  }
  
  # loop over methods
  results <- NULL
  for (method in methods) {
    if (is.null(results)) {
      results <- div_function(k, method, EM, ...)  # return results of the ret in a dataframe format if it is null
      colnames(results) <- "estimate"
    } else {
      next_results <- div_function(k, method, EM, ...)
      colnames(next_results) <- "estimate"
      results <- rbind.data.frame(results, next_results)  # if more than one method rbind all (not null)
    }
  }
  return(results)
}

calculateDiversity(k2, methods = c("simpson", "shannon", "invsimpson","chao1","chaobonge"))
