
#' Tree Dimension Test
#'
#' Computes the statistical significance for the presence of trajectory
#' in multivariate data.
#'
#' @param x matrix of input data. Rows as observations and columns as
#'  features.
#' @param perm number of simulations to compute null distribution
#'  parameters by maximum likelihood estimation.
#' @param MST the MST algorithm to be used in test. There are two
#'  options: "exact" MST and "boruvka" which is approximate but faster
#'  for large samples.
#' @param dim.reduction string parameter with value "pca" to perform
#'  dimensionality reduction or "none" to not perform dimensionality
#'  reduction before the test.
#'
#' @return A list with the following components:
#' \itemize{
#'  \item tdt_measure The tree dimension value for the given input data
#'  \item statistic  The S statistic calculated on the input data. S statistic is derived from tree dimension
#'  \item tdt_effect Effect size for tree dimension
#'  \item leaves Number of leaf/degree1 vertices in the MST  of the data
#'  \item diameter The tree diameter of MST, where each edge is of unit length
#'  \item p.value The pvalue for the S statistic. Pvalue measures presence of trajectory in input x.
#' }
#' @importFrom Rcpp evalCpp
#' @importFrom graphics plot legend
#' @importFrom stats dist plnorm rnorm cov prcomp
#' @importFrom mlpack emst
#' @import nFactors
#' @export
#' @useDynLib TreeDimensionTest

test.trajectory <- function(
  x, perm=100, MST=c("boruvka", "exact"),
  dim.reduction=c("pca", "none")
)
{

  #Error checking argument x
  if(missing(x)) stop("Argument x is required.")
  if(!is.matrix(x)) stop("Argument x should be a matrix.")
  if(perm < 0 ) stop("Argument perm should be greater thatn 0.")

  MST <- match.arg(MST)
  dim.reduction <- match.arg(dim.reduction)

  ###Apply pca if dim.reduction="pca"

  if(dim.reduction=="pca")
  {
    pc = prcomp(x)

    if(length(pc$sdev) >= 6) #nCng works with at least 6 variables
    {
      cmps = nCng(pc$sdev)
      x = pc$x[,c(1:cmps$nFactors)]
    }

    #if variables are less than 6, use all components
    else
    {
      x = pc$x
    }
  }


  #Getting null distribution
  dist = computeDists(x, perm, nrow(x), randomize.data, MST)


  #Getting observed statistic and tree dimension
  estimates = getStatistics(x, nrow(x), MST)

  #Fitting lognormal and calculating pvalue for s statistic
  fit = fitdist(dist, distr="lnorm", start = list(meanlog=1, sdlog =1))


  #Output result
  result = list()
  result$tdt_measure = estimates$td
  result$statistic = estimates$stat
  result$tdt_effect = estimates$effect
  result$leaves = estimates$leafs
  result$diameter = estimates$diameter
  result$p.value = plnorm(estimates$stat, meanlog = fit$estimate['meanlog'] , sdlog = fit$estimate['sdlog'], lower.tail = F)

  return(result)
}

#' Empirical Null Distribution of Tree Dimension Test
#'
#' Computes empirical null distribution of S statistic and parameters
#'  for lognormal approximation for input of size rows * columns using
#'  multivariate normal randomization

#' @param rows number of rows for data representing null case. Rows
#'  represent sample size.
#' @param cols number of columns for data representing null case.
#'  Columns represent variables.
#' @param perm number of simulations to compute null distribution.
#'  Default is 100.
#' @param MST name of MST to be used in computing distribution. There
#'  are two options; "exact" MST and "boruvka" which is faster for
#'  large samples

#' @return A list with the following components:
#' \itemize{
#'  \item dist       A vector with null distribution of s statistic
#'  \item meanlog    The meanlog parameter estimation for the lognormal distribution on empirical null distribution S.
#'  \item sdlog      The sdlog parameter estimation for lognormal distribution on empirical null distribution of S.
#' }
#' @import fitdistrplus
#' @export

empirical.distributions <- function(
  rows, cols, perm=100, MST=c("boruvka", "exact")
)
{
  if(missing(rows)) stop("Argument rows is required.")
  if(missing(cols)) stop("Argument cols is required.")
  if(perm < 0 ) stop("Argument perm should be greater thatn 0.")
  sample_size = rows
  x=randomize.data(rows,cols)
  distbns = computeDists(x, perm, sample_size, randomize.data, MST)


  #Fitting lognormal
  fit = fitdist(distbns, distr="lnorm", start = list(meanlog=1, sdlog =1))
  result = list()
  result$dist = distbns
  result$meanlog = fit$estimate['meanlog']
  result$sdlog = fit$estimate['sdlog']
  return(result)
}

#' Visualizing Euclidean Minimum Spanning Trees
#'
#' Plots an Euclidean minimum spanning tree from given input data. If
#'  labels for the samples in the data are provided, a legend is also
#'   plotted.

#' @param x input matrix, with rows as observations and columns as features
#' @param node.labels vector of labels for observations in x (vertices)
#' @param node.col vector of colors for the observations in x (vertices)
#' @param node.size numerical value to represent size of nodes in the plot
#' @param main title for the plot
#' @param legend.cord vector of the xy coordinates for the legend c(x,y)
#' @return result plots a minimum spanning tree for input data x
#' @import igraph
#' @import RColorBrewer
#' @export

plotTree <-  function(
  x, node.labels="", node.col="orange",
  node.size=5, main="MST plot",
  legend.cord=c(-1.2,1.1))
{

  if(missing(x)) stop("Argument x is required.")


  #Compute exact tree of x
  #mtree = make.tree(x)

  #transformation of data into igraph
  edges = get_edges(x, MST="boruvka")
  g=igraph::graph(edges, directed = FALSE)

  #coloring
  if(length(node.col)!=nrow(x))
  {

    #coloring
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    unique_labels = unique(node.labels)
    cols = vector(length=length(node.labels))
    j=1
    for(i in unique_labels)
    {
      cols[which(node.labels %in% i)]=j
      j=j+1
    }

    node.col = col_vector[cols]
  }

  plot(g, vertex.size=node.size ,vertex.color = node.col ,vertex.label="",main=main ,cex.main=0.03,layout=layout_with_kk, font.main=0.7)
  legend(legend.cord[1], legend.cord[2], legend = unique(node.labels) , col= unique(node.col), lwd=4, cex=0.6)

}



#Multivariate standard normal
randomize.data = function(sample, dim, iter, max_iter)
{
  #for reproducibility of null distribution
  set.seed(sample+dim+iter+max_iter)
  ran = rnorm(sample)

  for(i in c(2:dim))
  {
    ran = cbind(ran , rnorm(sample))
  }
  set.seed(NULL)
  return(ran)
}



#Get edges in correct format

# get_edges = function(m) {
#   edges = vector()
#
#   rows = sample(nrow(m))
#   cols = sample(ncol(m))
#   for(i in rows)
#   {
#     for(j in cols)
#     {
#       if(i<j && m[i,j])
#       {
#         edges = c(edges, i, j)}
#       }
#   }
#   return(edges)
#}


#' Tree Dimension Test Related Statistics
#'
#' Computes tree dimension measure, tree dimension test effect,
#' number leafs and tree diameter from MST of a given dataset

#' @param x matrix of input data. Rows as observations and columns as features
#' @param MST name of MST to be used in test. There are 2 options; "exact" MST and "boruvka" which is faster for large samples
#' @param dim.reduction string parameter with value "pca" to perform dimensionality reduction or "none" to not perform dimensionality reduction
#' @return A list with the following components:
#' \itemize{
#'  \item tdt_measure The tree dimension value for the given input data
#'  \item tdt_effect Effect size for tree dimension
#'  \item leaves Number of leaf/degree1 vertices in the MST  of the data
#'  \item diameter The tree diameter of MST, where each edge is of unit length
#'  }
#' @export


compute.stats <- function(
  x, MST=c("boruvka", "exact"),
  dim.reduction=c("pca", "none")
)
{
  MST <- match.arg(MST)
  dim.reduction <- match.arg(dim.reduction)

  if(dim.reduction=="pca")
  {
    pc = prcomp(x)

    if(length(pc$sdev) >= 6) #nCng works with at least 6 variables
    {
      cmps = nCng(pc$sdev)
      x = pc$x[,c(1:cmps$nFactors)]
    }

    #if variables are less than 6, use all components
    else
    {
      x = pc$x
    }
  }
  estimates = getStatistics(x, nrow(x), MST)
  result = list()
  result$tdt_measure = estimates$td
  result$tdt_effect = estimates$effect
  result$leaves = estimates$leafs
  result$diameter = estimates$diameter

  return (result)
}


#' Separability of Labeled Data Points
#'
#' Computes homogeneity of labeled observations with multiple label
#'  types.
#'
#' @param x input data matrix, with rows as observations and columns
#'  as features
#' @param labels a vector of labels for the observations. A label could
#'  be a type of the observation e.g cell type in single-cell data
#' @return A list with the following components:
#' \itemize{
#'  \item label_separability  A vector of separability scores for each
#'   of the label types. A high score denotes high separability
#'  \item overall_separability  Overall average separability score for
#'   all the labels
#' }
#'
#' @export
separability <- function(x,labels)
{
  unique_labels = unique(labels)
  label_purity = vector()
  edge_purity = vector()
  comb_purity = vector()
  tree = convert_to_tree(x)



  for(i in unique_labels)
  {

    #set of vertices for a unique label
    s = which(labels %in% i)

    #set starting index to 0
    s = s-1

    #Obtain minimum subtree cover of vertices in s
    res = minSubtreeCover(tree,s,labels)


    #Vertices in s
    a = length(s)

    #Label purity for each unique label
    label_purity[i] = a/res$subtree_nodes

  }


  return(list("label_separability"=label_purity, "overall_separability"= sum(label_purity)/length(unique_labels)))
}



make.tree <- function(x)
{
   dist = calcPairwiseDist(x)
   mtree = calculate_mst(dist)
  return(mtree)
}



