library(testthat)
library(LFQProcessing)

test_that("lfq des formatter works", 
          {
            des <- as.data.table(list(
              name = paste(c("OE0906_28331", "OE0906_28332", "OE0906_28333", 
                                  "OE0906_28334", "OE0906_28335", "OE0906_28336", "OE0906_28337"),".raw", sep=""),
              remove = c("False","False","False","False","False","False", "True"),
              bioRep = seq(7),
              techRep = seq(7),
              mqExperiment = c("018DBS2","013DBS1","021DBS2","018DBS1","024DBS2","012DBS2B","to_be_removed"),
              experiment = c("A","A","A","B","B","B",""))
            )

            des_list <- LFQProcessing::experiment_design_formatter(des)
            des <- des_list[[1]]

            expect_equal(des$file_name,  c("OE0906_28331", "OE0906_28332", "OE0906_28333", 
                                           "OE0906_28334", "OE0906_28335", "OE0906_28336"))
            expect_equal(nrow(des),6 )
            expect_equal(unique(des$remove), "FALSE")
            expect_false("to_be_removed" %in% des$mqExperiment)
            
          }
          
)
