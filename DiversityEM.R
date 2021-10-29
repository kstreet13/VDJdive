#EMDiversity
createCountsMatrixEm<-function(x) {
  df<-as.data.frame(colSums(x$clono))
  return(df)
}

k2<-createCountsMatrixEm(sce)


createCountsMatrixEm(sce)

calculateDiversity <- function(k, methods, EM = FALSE , ...) { #any parm but k or methods is ... for e.g. cutoff for breakaway which will pass down
  # function to calculate diversity for a given method 
  div_function <- function(k, method,EM, ...) {
    if (method %in% c("simpson", "shannon", "invsimpson")) {
      # simpson, shannon and invsimpson calculated using diversity()
      ret <- data.frame(t(diversity(t(k), ..., index=method)))
      rownames(ret) <- paste0(method, ".estimate") #paste .est after running fun for rownames
    } else if ("chao1" == method & !EM) {
      # use fossil::chao1() for chao1
      ret <- data.frame(t(apply(k, 2, fossil::chao1, ...))) #apply chao1 fun for every column
      rownames(ret) <- paste0(method, ".estimate")
    } else if ("chao1" == method  & EM) {
      # use fossil::chao1() for chao1
      ret <- round(k*100,0)
      ret1 <- data.frame(t(apply(ret, 2, fossil::chao1, ...))) #apply chao1 fun for every column
      rownames(ret1) <- paste0(method, ".estimate")
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
      ret1 <- data.frame(apply(ret, 2, function(x, ...) {
        bw <- breakaway(x, ...)
        c(estimate=bw$estimate, error=bw$error)
      }))
      rownames(ret1) <- paste0(method, ".", rownames(ret1))
    }
    else {
      stop ("invalid diversity method")
    } #error message we can change later
    return(ret1)
  }
  
  # loop over methods
  results <- NULL
  for (method in methods) {
    if (is.null(results)) {
      results <- div_function(k, method, EM, ...)#return results of the ret in a dataframe format if it is null
    } else {
      results <- rbind.data.frame(results, div_function(k, method, ...)) #if more than one method rbind all (not null)
    }
  }
  return(results)
}




