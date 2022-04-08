library(LFQProcessing)
library(testthat)
library(tools)


test_that("Returns lfq intensity when all intensity columns are present", 
          {
            
            des = as.data.table(list(
              file_name = c("A","B","C"),
              fraction = c(0,0,0),
              experiment = c("bdcya","bdcya","bdcya"),
              remove = c(FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,2,3)
            ))
            
            dt = as.data.table(list(
              fasta_headers = c("blablabla","blablabla2"),
              `intensity sample1rep1` = c(1,1),
              `intensity sample1rep2` = c(1,1),
              `intensity sample1rep3` = c(1,1),
              `lfq intensity sample1rep1` = c(1,1),
              `lfq intensity sample1rep2` = c(1,1),
              `lfq intensity sample1rep3` = c(1,1)
            )
            )

            lfq_intensity_cols = get_lfq_intensity_columns(dt, des)
            expect_equal(lfq_intensity_cols,
                         c("intensity sample1rep1", "intensity sample1rep2", "intensity sample1rep3",
                           "lfq intensity sample1rep1", "lfq intensity sample1rep2", "lfq intensity sample1rep3"))
          })

test_that("Returns intensity columns when all lfq intensity columns not present",
          {
          des = as.data.table(list(
            file_name = c("A","B","C"),
            fraction = c(0,0,0),
            experiment = c("bdcya","bdcya","bdcya"),
            remove = c(FALSE,FALSE,FALSE),
            mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
            Replicate = c(1,2,3)
          ))
          
          dt = as.data.table(list(
            fasta_headers = c("blablabla","blablabla2"),
            `intensity sample1rep1` = c(1,1),
            `intensity sample1rep2` = c(1,1),
            `intensity sample1rep3` = c(1,1)
          )
          )

          lfq_intensity_cols = get_lfq_intensity_columns(dt, des)
          expect_equal(lfq_intensity_cols,
                       c("intensity sample1rep1", "intensity sample1rep2", "intensity sample1rep3"))
          }
          )

test_that("Returns intensity columns when some lfq intensity columns not present",
          {
          des = as.data.table(list(
            file_name = c("A","B","C"),
            fraction = c(0,0,0),
            experiment = c("bdcya","bdcya","bdcya"),
            remove = c(FALSE,FALSE,FALSE),
            mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
            Replicate = c(1,2,3)
          ))
          
          dt = as.data.table(list(
            fasta_headers = c("blablabla","blablabla2"),
            `intensity sample1rep1` = c(1,1),
            `intensity sample1rep2` = c(1,1),
            `intensity sample1rep3` = c(1,1),
            `lfq intensity sample1rep2` = c(1,1),
            `lfq intensity sample1rep3` = c(1,1)
          )
          )
          lfq_intensity_cols = get_lfq_intensity_columns(dt, des)
          expect_equal(lfq_intensity_cols,
                       c("intensity sample1rep1", "intensity sample1rep2", "intensity sample1rep3"))
          }
          )

test_that("Error when missing some intensity columns are required", 
          {
          des = as.data.table(list(
            file_name = c("A","B","C"),
            fraction = c(0,0,0),
            experiment = c("bdcya","bdcya","bdcya"),
            remove = c(FALSE,FALSE,FALSE),
            mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
            Replicate = c(1,2,3)
          ))
          
          dt = as.data.table(list(
            fasta_headers = c("blablabla","blablabla2"),
            `intensity sample1rep1` = c(1,1),
            `intensity sample1rep3` = c(1,1),
            `lfq intensity sample1rep2` = c(1,1),
            `lfq intensity sample1rep3` = c(1,1)
          )
          )

          lfq_intensity_cols = get_lfq_intensity_columns(dt, des)
          expect_equal(lfq_intensity_cols, NA)
          }
          )

test_that("LFQ Fractions gets correct columns (duplicate in des)", 
          {
            
            des = as.data.table(list(
              file_name = c("A","B","C"),
              fraction = c(0,0,0),
              experiment = c("bdcya","bdcya","bdcya"),
              remove = c(FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep1","sample1rep3"),
              Replicate = c(1,1,3)
            ))
            
            dt = as.data.table(list(
              fasta_headers = c("blablabla","blablabla2"),
              `intensity sample1rep1` = c(1,1),
              `intensity sample1rep3` = c(1,1),
              `lfq intensity sample1rep1` = c(1,1),
              `lfq intensity sample1rep3` = c(1,1)
            )
            )
            
            lfq_intensity_cols = get_lfq_intensity_columns(dt, des)
            expect_equal(lfq_intensity_cols,
                         c("intensity sample1rep1", "intensity sample1rep3",
                           "lfq intensity sample1rep1", "lfq intensity sample1rep3"))
          })

test_that("LFQ with samples removed works", 
          {
            
            des = as.data.table(list(
              file_name = c("A","B","C"),
              fraction = c(0,0,0),
              experiment = c("bdcya","bdcya","bdcya"),
              remove = c(FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2"),
              Replicate = c(1,2,3)
            ))
            
            dt = as.data.table(list(
              fasta_headers = c("blablabla","blablabla2"),
              `intensity sample1rep1` = c(1,1),
              `intensity sample1rep2` = c(1,1),
              `intensity sample1rep3` = c(1,1),
              `lfq intensity sample1rep1` = c(1,1),
              `lfq intensity sample1rep2` = c(1,1),
              `lfq intensity sample1rep3` = c(1,1)
            )
            )
            
            lfq_intensity_cols = get_lfq_intensity_columns(dt, des)
            expect_equal(lfq_intensity_cols,
                         c("intensity sample1rep1", "intensity sample1rep2",
                           "lfq intensity sample1rep1", "lfq intensity sample1rep2"))
          })

