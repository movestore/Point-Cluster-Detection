library("move2")

test_data <- test_data("input3_move2.rds") # file must be move2!

test_that("dur_unit mins executes", {
  actual <- rFunction(data = test_data, dur_unit = "mins")
  expected_count <- 821
  expect_equal(nrow(actual), expected_count)
})


test_that("gap_unit mins executes", {
  actual <- rFunction(data = test_data, gap_unit = "mins")
  expected_count <- 0
  expect_equal(nrow(actual), expected_count)
})


test_that("dur_unit non mins executes", {
  actual <- rFunction(data = test_data, dur_unit = "days")
  expected_count <- 359
  expect_equal(nrow(actual), expected_count)
})


test_that("gap_unit non mins executes", {
  actual <- rFunction(data = test_data, gap_unit = "days")
  expected_count <- 359
  expect_equal(nrow(actual), expected_count)
})
