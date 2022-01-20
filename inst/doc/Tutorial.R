## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE,message=FALSE,
  collapse = TRUE,
  comment = "#>"
  
)

## ----setup--------------------------------------------------------------------
library(TreeDimensionTest)

## -----------------------------------------------------------------------------
data_path=system.file('extdata', package = 'TreeDimensionTest')
input = readRDS(paste0(data_path,"/bifurcating_7.rds"))

#'input' is matrix of single-cell RNA-seq gene expression data with rows as cells and columns as genes
input = input$expression

#PCA and plot to visualize data
pc=prcomp(input)
plot(pc$x[,c(1:2)], main="PCA plot", xlab="PC1", ylab="PC2")

#Running the test.trajectory function on matrix "input".

# dim.reduction is set to "pca"; meaning dimensionality reduction will be performed first #using principal component analysis.Number of pca components are selected using Scree test. Set dim.reduction to "none" if you don't wish to perform #dimensionality reduction.

#MST is set to "exact"; the exact MST is used. The alternative is the approximate and fast
#DualTreeBoruvka MST. Set MST to "boruvka" to use the approximate MST.
res = test.trajectory(input, dim.reduction = "pca", MST="exact")

#List containing Tree Dimension Test measure (tdt_measure), Tree Dimension Test #effect(tdt_effect), statistic, p.value, number vertices that are leaves and diameter of tree.
#p.value is significant and tdt_effect is strong, depicting presence of trajectory.
res

## -----------------------------------------------------------------------------
mat = cbind(rnorm(1000), rnorm(1000))
res = test.trajectory(mat, dim.reduction = "none")
plot(mat)
res

## -----------------------------------------------------------------------------
mat = cbind(rnorm(1000), rnorm(1000))

res = compute.stats(mat, MST="boruvka", dim.reduction = "none")
res

## -----------------------------------------------------------------------------
#Random data
mat = cbind(rnorm(200), rnorm(200))

#Labels for the samples in the data
labels = c(rep("L1", 93), rep("L2",78), rep("L3",29))

#Color vector of samples, each unique color correspods with unique label
cols = c(rep("blue",93), rep("green",78), rep("red",29))

#Plots an MST of the data, with samples of the same label highlighted by same color
plotTree(mat,labels, node.size = 12, node.col = cols,main = "Low seperability", legend.cord=c(-2.1,0.9))

#Compute separability of samples in mat
res = separability(mat, labels)

#List containing separability values for each label and, overall separability on the data. #Overall separability is relatively low, implying samples with same labels are mixed.
res

## -----------------------------------------------------------------------------
#An example data where labels of the same type are close together , resulting in high separability value.
mat = rbind(cbind(rnorm(93,mean=20), rnorm(93, mean=20)), cbind(rnorm(78,mean=5),rnorm(78,mean=5)), cbind(rnorm(29, mean=50), rnorm(29, mean=50)))
labels = c(rep("L1", 93), rep("L2",78), rep("L3",29))

#Color vector of samples, each unique color correspods with unique label
cols = c(rep("blue",93), rep("green",78), rep("red",29))

plotTree(mat,labels, node.size=12, node.col = cols, main = "High seperability", legend.cord=c(-1.9,0.9))
res = separability(mat, labels)

#Overall separability is 1, implying  samples of different labels are perfectly separated.
res

## -----------------------------------------------------------------------------
#Loading calcium signaling pathway data from Mouse
#This is Mouse development RNA-seq data spanned by the geneset fo calcium signaling pathway
#Rows are genes and columns are samples.
data_path=system.file('extdata', package = 'TreeDimensionTest')
load(file=paste0(data_path,"/calcium_pathway_data.rdata"))

#loading color vector of samples by label type; mouse_cols
load(file=paste0(data_path,"/mouse_cols.rdata"))

#Labels of the samples are the column names of the data, which are names of tissue types.
labels = colnames(calcium_pathway_data)

plotTree(t(calcium_pathway_data), labels, node.col=mouse_cols,node.size=12, main = "Calcium Signaling pathway", legend.cord=c(-1.9,-1.3))


res = separability(t(calcium_pathway_data), labels)

#Separabiltiy for each tissue type as well as the overall separability. High separability depicts high tissue specificity.
res



## -----------------------------------------------------------------------------
#Loading ribosome pathway data from Mouse
data_path=system.file('extdata', package = 'TreeDimensionTest')
load(file=paste0(data_path,"/ribosome_pathway_data.rdata"))

#loading color vector of samples by label type; mouse_cols
load(file=paste0(data_path,"/mouse_cols.rdata"))

# ribsome_pathway_data is RNA-seq data with rows as genes and columns as samples
labels = colnames(ribosome_pathway_data)
plotTree(t(ribosome_pathway_data), labels, node.col= mouse_cols, node.size=12, main = "Ribosome pathway", legend.cord=c(-1.9,-1.3))
res = separability(t(ribosome_pathway_data), labels)
res


