#TDT on linear tree
x = rnorm(100)
y = x+1
linear = cbind(x,y)

test_that("Accuracy of tree dimension test (test.trajectory) statistics on linear MST", {
  expect_equal(test.trajectory(linear, MST="exact", dim.reduction = "none")$tdt_measure, 1)
  expect_equal(test.trajectory(linear, MST="exact", dim.reduction = "none")$tdt_effect, 1)
  expect_equal(test.trajectory(linear, MST="exact", dim.reduction = "none")$leaves, 2)
  expect_equal(test.trajectory(linear, MST="exact", dim.reduction = "none")$diameter, 99)
})


test_that("Accuracy of tree dimension test (compute.stats) statistics on linear MST", {
  expect_equal(compute.stats(linear, MST="exact", dim.reduction = "none")$tdt_measure, 1)
  expect_equal(compute.stats(linear, MST="exact", dim.reduction = "none")$tdt_effect, 1)
  expect_equal(compute.stats(linear, MST="exact", dim.reduction = "none")$leaves, 2)
  expect_equal(compute.stats(linear, MST="exact", dim.reduction = "none")$diameter, 99)
})


# Tdt on star tree
x = c(0,0,0,1,-1)
y = c(0,1,-1,0,0)
star = cbind(x,y)


test_that("Accuracy of tree dimension test (test.trajectory) statistics on star MST", {
  expect_equal(compute.stats(star, MST="exact", dim.reduction = "none")$tdt_measure, 2)
  expect_equal(compute.stats(star, MST="exact", dim.reduction = "none")$tdt_effect, 0)
  expect_equal(compute.stats(star, MST="exact", dim.reduction = "none")$leaves, 4)
  expect_equal(compute.stats(star, MST="exact", dim.reduction = "none")$diameter, 2)
})



##Dimension 2 case; TDT on dimension 2 tree
x = seq(from=0, to = 10, by=1)
y = rep(5,length(x))
mat = cbind(x,y)

mat = rbind(mat, cbind(rep(5,length(x)),x))
mat = mat[-6,]

test_that("Accuracy of tree dimension test (test.trajectory) statistics on linear MST", {
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$tdt_measure, 2)
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$tdt_effect, 0.699, tolerance=0.001)
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$leaves, 4)
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$diameter, 10)
})


test_that("Accuracy of tree dimension test (test.trajectory) statistics on linear MST", {
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$tdt_measure, 2)
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$tdt_effect, 0.699, tolerance=0.001)
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$leaves, 4)
  expect_equal(compute.stats(mat, MST="boruvka", dim.reduction = "none")$diameter, 10)
})
