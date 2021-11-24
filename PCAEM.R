#EMDiversity

library(breakaway)
library(fossil)
library(vegan)

FileToy<-readRDS("~/Desktop/PostDoColabProj/toyTCRdata.rds")
quant_Toy<-EMquant(FileToy)


#Function to create a matrix compatibe for diversity function
createCountsMatrixEm<-function(ranged_sum_exp) {
  col_data <- colData(ranged_sum_exp) #EM_alg data matrix
  by_sample <- split(col_data, col_data$sample) # split matrix by samples 
  col_sums <- lapply(by_sample, function(df) colSums(df$clono)) #perform colsums for colontypes per samples
  final_df <- as.data.frame(col_sums) #make list to dataframe
  return(final_df)
}

k2<-createCountsMatrixEm(quant_Toy)

#PCA function compatible createCountMatrixEM output

PcaTcrSeq<-function(k, application ="Sample") {
      if(application == "TCR"){
        PCA_ret <- prcomp(t(k))
      } else if (application == "Sample") {
        PCA_ret <- prcomp(k) 
      } else{
        stop("Invalid value provided for application, the value should be either `TCR` or `Sample`")
      }
    return(PCA_ret)
}

pcaobj<-PcaTcrSeq(k2,"Sample")   

library(ggplot2)
library(RColorBrewer)

#PlotPCA rotation (samples)
PlotPCA<-function(pcaobj) {
        plotDATA1 <- data.frame(pcaobj$rotation,
                        headers=rownames(pcaobj$rotation),
                        color=apply(pcaobj$rotation, 1, function(x) sum(x)))
        qplotObj<-qplot(PC1,PC2,data=plotDATA1,label=rownames(plotDATA1), color=color)+
                  geom_point(size=2) + 
                  ggtitle("PCA of data")+
                  theme(plot.title = element_text(size = 10, face = "bold"))
        qplotObj<-qplotObj + scale_color_gradientn(colours = rainbow(5))
        return(qplotObj)
}

PlotPCA(pcaobj)

