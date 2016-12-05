context("bootstrap arguments")
test_that("check_matrix with ok vectors", {
  expect_true(check_data(matrix(1:9,3,3)))
})
test_that("check_matrix fails with invalid argument",{
  expect_error(check_data(list(x=c(1:3),y=c(2:6))))
  expect_error(check_data(matrix(c("a","b","c","d"),2,2)))
})
test_that("check_sets with ok sets", {
  expect_true(check_sets(list(1:3,4:5)))
  expect_true(check_sets(list("a","b")))
})
test_that("check_sets with invalid argument", {
  expect_error(check_sets(list(1:3,"c")))
})
test_that("check_ncomps with ok sets", {
  expect_true(check_ncomps(data=matrix(1:9,3,3),2))
})
test_that("check_ncomps with invalid argument", {
  expect_error(check_ncomps(data=matrix(1:9,3,3),"c"))
  expect_error(check_ncomps(data=matrix(1:9,3,3),c(1,2)))
  expect_error(check_ncomps(data=matrix(1:9,3,3),4))
})
test_that("check_numeric with ok vectors", {
  expect_true(check_numeric(matrix(1:9,3,3)))
})
test_that("check_numeric fails with invalid argument",{
  expect_error(check_numeric(matrix(c("a","b","c","d"),2,2)))
})


