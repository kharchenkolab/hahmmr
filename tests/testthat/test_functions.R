
library(hahmmr)

test_that("dpoilog works", {
  expect_equal(hahmmr:::dpoilog(c(1,11),c(1,1),c(1,1)), c(0.175733342664327, 0.0150105250670325))
})

test_that("logSumExp() works", {

  b = hahmmr:::logSumExp(c(1.2, 3.4))
  expect_equal(b, 3.5050833)  
  d = hahmmr:::logSumExp(c(1.2, 3.4, -5.6, -7.8))
  expect_equal(d, 3.5052067) 
   
})

test_that("Check that likelihood_allele() works as expected", {

  LL = hahmmr:::likelihood_allele(pre_likelihood_hmm)
  expect_equal(as.integer(LL), -736)

})

test_that("Check that forward_backward() works as expected", {

  p_major = hahmmr:::forward_back_allele(pre_likelihood_hmm)[,1]
  expect_equal(is.vector(p_major), TRUE)
  expect_equal(length(p_major), 1042)
  expect_equal(round(p_major[1], 3), 0.963)
  expect_equal(round(p_major[2], 3), 0)
  expect_equal(round(p_major[3], 3), 0)
  expect_equal(round(p_major[10], 3), 0.745)

})

test_that("Check that viterbi_allele() works as expected", {

  states = hahmmr:::viterbi_allele(pre_likelihood_hmm)
  expect_equal(sum(states == 1), 440)
  expect_equal(sum(states == 2), 602)
  expect_equal(states[1], 1)
  expect_equal(states[2], 2)
  expect_equal(states[3], 2)
  expect_equal(states[1024], 1)
  
})
