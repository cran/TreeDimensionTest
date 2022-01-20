# completely separated
x = seq(from=1, to=20, by=1)
y = 0.1 * x * x
mat = cbind(x,y)
labels = c(rep("red",12), rep("blue",8))

test_that("Separability on linear completely separated MST", {
  res = separability(mat, labels=labels)
  expect_equal(unname(res$label_separability['red']), 1)
  expect_equal(unname(res$label_separability['blue']), 1)
  expect_equal(res$overall_separability, 1)
})



# Mixed, no adjacent samples with same label

x = seq(from=1, to=20, by=1)
y = 0.1 * x * x
mat = cbind(x,y)
labels = rep(c("red","blue"),10)

test_that("Separability on linear completely mixed MST", {
  res = separability(mat, labels=labels)
  expect_equal(unname(res$label_separability['red']), 10/19)
  expect_equal(unname(res$label_separability['blue']), 10/19)
  expect_equal(res$overall_separability, 10/19)
})
