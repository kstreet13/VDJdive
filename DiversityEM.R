library(sparseMatrixStats)
library(dplyr)
library(breakaway)
require(SingleCellExperiment)
sce <- readRDS('~/Desktop/PostDoColabProj/toyTCRdata.rds')
source('EMquant.R')
sce <- EMquant(sce)
# the resulting matrix:
sce$clono 
table(sce$clono@i)
table(sce$clono@p)
table(sce$clono@Dim)
table(sce$clono@x)

#Split By Sample and run diveristy

#thing 1 split samples before you run the diversity measures


EM_Div_input<-as.data.frame(colSums(sce$clono))
kk12<-as.data.frame(sample(0.0:50.1, 29322, replace = TRUE)) 
EM_Div_input2<-cbind(EM_Div_input,kk12)
names(EM_Div_input2)<-c("ini","non")
calculate_diversity(EM_Div_input2,"shannon")
calculate_diversity(EM_Div_input2,"simpson")
calculate_diversity(EM_Div_input2,"invsimpson")
#for Chao1 bundge int is needed
EM_Div_input3<-EM_Div_input2
EM_Div_input3<-round(EM_Div_input3*10,0)
calculate_diversity(EM_Div_input3,"chao1")
calculate_diversity(EM_Div_input3,"chaobonge")

#Diversity for EM
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




