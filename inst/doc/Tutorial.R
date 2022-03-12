## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE,
  collapse = TRUE, cache = FALSE, comment = "#>", out.width = "65%")

## -----------------------------------------------------------------------------
require(TreeDimensionTest)

data_path=system.file('extdata', 'bifurcating_7.rds', package = 'TreeDimensionTest')

input = readRDS(data_path)

## -----------------------------------------------------------------------------
res = test.trajectory(input)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(unlist(res[-7], recursive = FALSE), col.names = "Test results contained in `res`")

## -----------------------------------------------------------------------------
plot(res, node.size=12, node.col="mediumseagreen", 
     main = paste0("Trajectory presence\np-value = ",
                   format.pval(res$p.value, digit=2)))
pc=prcomp(input)
plot(pc$x[,c(1:2)], xlab="PC1", ylab="PC2", cex=1,
     sub="Principal components", col="mediumseagreen",
     main = paste0("Trajectory presence\np-value = ",
                   format.pval(res$p.value, digit=2), 
                   "\n(n = ", nrow(input), ")")
    )

## -----------------------------------------------------------------------------
res = test.trajectory(input, dim.reduction = "pca", MST = "exact")

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(unlist(res[-7], recursive = FALSE), col.names = "Test results contained in `res`")

## -----------------------------------------------------------------------------
n = 100
mat = cbind(rnorm(n), rnorm(n))
res = test.trajectory(mat, dim.reduction = "none")

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(unlist(res[-7]), col.names = "Test results contained in `res`")

## -----------------------------------------------------------------------------
plot(res, node.size=12, 
     main = paste0("Trajectory presence\np-value = ",
                   format.pval(res$p.value, digit=2)))
plot(mat, col="orange", pch=2, cex=0.5, xlab="Dim 1", ylab="Dim 2",
     main=paste0("Trajectory presence\np-value = ", 
                format.pval(res$p.value, digit=2), 
                "\n(n = ", n, ")")
     )

## -----------------------------------------------------------------------------
res = compute.stats(mat, MST="boruvka", dim.reduction = "none")

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(unlist(res[-5]), col.names = "Results contained in `res`")

## -----------------------------------------------------------------------------
#Random data
mat = cbind(rnorm(200), rnorm(200))

#Labels for the samples in the data
labels = c(rep("L1", 93), rep("L2",78), rep("L3",29))

#Color vector of samples, each unique color correspods with unique label
cols = c(rep("orange", 93), rep("mediumseagreen", 78), rep("purple", 29))

#Compute separability of samples in mat
res = separability(mat, labels)

## -----------------------------------------------------------------------------
knitr::kable(res$label_separability, col.names = "Label separability")

## -----------------------------------------------------------------------------
#Plots an MST of the data, with samples of the same label highlighted by same color
# plotTree(mat, labels, node.size = 12, node.col = cols, 
#   main = paste("Low seperability", format(res$overall_separability, digits = 3)), 
#   legend.cord=c(-2.1,0.9)
# )

plot(res,node.size = 12, node.col = cols, 
  main = paste("Low seperability", format(res$overall_separability, digits = 3)), legend.cord=c(-2.1,0.9))

## -----------------------------------------------------------------------------
mat = rbind(
  cbind(rnorm(93,mean=20), rnorm(93, mean=20)),
  cbind(rnorm(78,mean=5), rnorm(78,mean=5)),
  cbind(rnorm(29, mean=50), rnorm(29, mean=50)))
labels = c(rep("L1", 93), rep("L2",78), rep("L3",29))
res = separability(mat, labels)


## -----------------------------------------------------------------------------
knitr::kable(res$label_separability, col.names = " Label separability")

## -----------------------------------------------------------------------------
#Color vector of samples corresponding to labels
cols = c(rep("orange", 93), rep("mediumseagreen", 78), rep("purple", 29))

# plotTree(
#   mat,labels, node.size=12, node.col = cols, 
#   main = paste("High seperability", format(res$overall_separability, digits = 3)), 
#   legend.cord=c(-1.9,0.9))


plot(res,node.size = 12, node.col = cols, 
  main = paste("High seperability", format(res$overall_separability, digits = 3)), legend.cord=c(-2.1,0.9))

## -----------------------------------------------------------------------------
#Loading calcium signaling pathway data from Mouse
#This is Mouse development RNA-seq data spanned by the geneset fo calcium signaling pathway
#Rows are genes and columns are samples.

file.rdata=system.file('extdata', 'calcium_pathway_data.rdata', package = 'TreeDimensionTest')
load(file=file.rdata)

#loading color vector of samples by label type; mouse_cols
file.rdata=system.file('extdata', 'mouse_cols.rdata', package = 'TreeDimensionTest')
load(file=file.rdata)

#Labels of the samples are the column names of the data, which are names of tissue types.
labels = colnames(calcium_pathway_data)

res = separability(t(calcium_pathway_data), labels)

#Separabiltiy for each tissue type as well as the overall separability. High separability depicts high tissue specificity.


# plotTree(
#   t(calcium_pathway_data), labels, node.col=mouse_cols,node.size=12, 
#   main = paste("Calcium signaling pathway\nTissue specificity", 
#                format(res$overall_separability, digits = 3)), 
#   legend.cord=c(-1.9,-1.3))

plot(res,node.size = 12, node.col=mouse_cols,
 main = paste("Calcium signaling pathway\nTissue specificity", 
               format(res$overall_separability, digits = 3)), legend.cord=c(-1.9,-1.3))

## -----------------------------------------------------------------------------
knitr::kable(res$label_separability, col.names = "Calcium signaling pathway tissue specificity")

## -----------------------------------------------------------------------------
#Loading ribosome pathway data from Mouse
file.rdata=system.file('extdata', 'ribosome_pathway_data.rdata', package = 'TreeDimensionTest')
load(file=file.rdata)

#loading color vector of samples by label type; mouse_cols
file.rdata=system.file('extdata', 'mouse_cols.rdata', package = 'TreeDimensionTest')
load(file=file.rdata)

# ribsome_pathway_data is RNA-seq data with rows as genes and columns as samples
labels = colnames(ribosome_pathway_data)
res = separability(t(ribosome_pathway_data), labels)
plot(
  res, node.col= mouse_cols,
  node.size=12, 
  main = paste("Ribosome pathway\nTissue specificity",
               format(res$overall_separability, digits = 3)), 
  legend.cord=c(-1.9,-1.3))

## -----------------------------------------------------------------------------
knitr::kable(res$label_separability, col.names = "Ribosome pathway tissue specificity")

