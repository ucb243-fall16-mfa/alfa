context("plot arguments")
test_that("check_type with ok vectors", {
  expect_true(check_type(1))
})
test_that("check_type fails with invalid argument",{
  expect_error(check_type("a"))
  expect_error(check_type(c(1,2)))
  expect_error(check_type(5))             
})
test_that("check_scale with ok vectors", {
  expect_true(check_scale(1,2))
})
test_that("check_scale fails with invalid argument",{
  expect_error(check_scale(1,"c"))
  expect_error(check_scale(1,c(1:2)))
})