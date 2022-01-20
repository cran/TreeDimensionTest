NEWS
====
2022-01-17
1. Replaced the "ComputeMST" function from emstreeR package with "emst" function from mlpack package to compute fast MST.

2. Modified "getStatistics" function to accommodate output of "emst" function which is now a list. Output of "ComputeMST" was a dataframe.

3.Removed 'randomize' function
4.Removed 'DFtoNM' function

2021-12-12
1. Created inst/extdata to store data used by vignette. Added data files bifurcating_7.rds,       calcium_pathway_data.rdata and ribosome_pathway_data.rdata.

2. Created testcases using testthat package. Created test-testrajectory.R file, which contains    test cases for functions 'test.trajectory' and 'compute.stats'. Also created
   test-separability.R file which contains test cases for function 'separability'.

3. Modified function 'plotTree' to accept node colors as argument

4. Modified vignette to improve visibility of separability plots 

2021-12-11

1. TreeDimensionTest version 0.0.1 created from version 1.0.
2. Replaced depreciated random_shuffle() by shuffle() in tree_dim_src.cpp.

2021-12-01

1. TreeDimensionTest version 1.0 created.
