---
title: NEWS
---

## Version 0.0.2

2022-03-11
  
  1. Created package version 0.0.2 from 0.0.1.1.

## Version 0.0.1.1 (not publicly released)

2022-03-10

  1. Removed the `method` parameter in `separability()` function. Separability is computed using the formula presented in the paper
  2. Modified `minSubtreeCover()` function in tree_dim_src.cpp to represent a tree by an adjacency list.
  3. Included code to check for required parameters in `separability()` function.


2022-02-10

  1. Added a `method` parameter to the `separability()` function. A value of `"original"` returns tissue specificity as used in the original paper. A value of `"normalized"` returns the normalized tissue specificity.
  2. Modified `test.trajectory()` and `compute.stats()` to return size of original dimension and number of PCA components.
  3. Modified `test.trajectory()` to return a class of type `"treedim"` and to return an additional component `"mst"` in the list.
  4. Modified `separability()` to return a class of type `"treedim"` and to return an additional component `"mst"` in the list.
  5. Created a function `plot.treedim()` that plots an MST for an object of type `"treedim"`.


2022-01-20
  
  1. Created package version 0.0.1.1 from 0.0.1.
  2. License changed from GPL (>= 2) to LGPL (>= 3).
  3. Updated CITATION.
  4. Included "@importFrom Rdpack reprompt" to avoid NOTE to 0.0.1 on CRAN.
  
## Version 0.0.1

2022-01-17

  1. Replaced the "ComputeMST" function from emstreeR package with "emst" function from mlpack package to compute fast MST.

  2. Modified "getStatistics" function to accommodate output of "emst" function which is now a list. Output of "ComputeMST" was a dataframe.

  3. Removed 'randomize' function

  4. Removed 'DFtoNM' function

2021-12-12

  1. Created inst/extdata to store data used by vignette. Added data files bifurcating_7.rds, calcium_pathway_data.rdata and ribosome_pathway_data.rdata.

  2. Created testcases using testthat package. Created test-testrajectory.R file, which contains test cases for functions 'test.trajectory' and 'compute.stats'. Also created test-separability.R file which contains test cases for function 'separability'.

  3. Modified function 'plotTree' to accept node colors as argument

  4. Modified vignette to improve visibility of separability plots 

2021-12-11

  1. TreeDimensionTest version 0.0.1 created from version 1.0.
  2. Replaced depreciated random_shuffle() by shuffle() in tree_dim_src.cpp.

## Version 1.0

2021-12-01

  1. TreeDimensionTest version 1.0 created.
