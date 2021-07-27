library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)


acceptance_test<- function(current, expected,
                           tolerance = 10**-3){
  
  test_that("approximately equal", {
    
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
}

data("example_tmt_PXD019880")

current_tmt <- tmt_transformer(protein_groups,
                               des,
                               "./tmp",
                               imputeStDev=0.3,
                               imputePosition=1.8)

#plot(current_tmt[[1]]$`logFC Cerebellum - Cerebrum`, -log10(current_tmt[[1]]$`P.Value Cerebellum - Cerebrum`))


load("../data/expected_tmt_PXD019880.rda")

acceptance_test(current_tmt[[1]], expected_tmt[[1]])