context("Lg arguments")
test_that("check_numeric with ok vectors", {
  expect_true(check_numeric(matrix(1:9,3,3)))
})
test_that("check_numeric fails with invalid argument",{
  expect_error(check_numeric(matrix(c("a","b","c","d"),2,2)))
})
