library(LFQProcessing)
library(testthat)
library(tools)

# replicate	reporter_channel	condition
# 1	0	A
# 1	1	B
# 2	2	A
# 2	3	B
# 3	4	A


test_that("Returns tmt intensity when there is file name", 
          {
            
            des = as.data.table(list(
              experiment = c("EV", "EV", "EV", "EV", "EV"),
              condition = c("A","B","A", "B","A"),
              fraction = c(0,0,0,0,0,0,0),
              experiment = c("","","","",""),
              remove = c(FALSE,FALSE,FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,1,2,2,3),
              reporter_channel = c(1,2,3,4,5)
            ))
            
            dt = as.data.table(list(
              fasta_headers = c("blablabla","blablabla2"),
              `reporter intensity corrected 1 EV` = c(1,1),
              `reporter intensity corrected 2 EV` = c(1,1),
              `reporter intensity corrected 3 EV` = c(1,1),
              `reporter intensity corrected 4 EV` = c(1,1),
              `reporter intensity corrected 5 EV` = c(1,1)
            )
            )
            
            tmt_intensity_cols = get_tmt_intensity_columns(dt, des)
            expect_equal(tmt_intensity_cols,
                         c("reporter intensity corrected 1 EV",
                           "reporter intensity corrected 2 EV",
                           "reporter intensity corrected 3 EV",
                           "reporter intensity corrected 4 EV", 
                           "reporter intensity corrected 5 EV"))
          })

test_that("Returns tmt intensity when there is file name extra space", 
          {
            
            des = as.data.table(list(
              experiment = c("EV", "EV", "EV", "EV", "EV"),
              condition = c("A","B","A", "B","A"),
              fraction = c(0,0,0,0,0,0,0),
              experiment = c("","","","",""),
              remove = c(FALSE,FALSE,FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,1,2,2,3),
              reporter_channel = c(1,2,3,4,5)
            ))
            
            dt = as.data.table(list(
              fasta_headers = c("blablabla","blablabla2"),
              `reporter intensity corrected 1  EV` = c(1,1),
              `reporter intensity corrected 2  EV` = c(1,1),
              `reporter intensity corrected 3  EV` = c(1,1),
              `reporter intensity corrected 4  EV` = c(1,1),
              `reporter intensity corrected 5  EV` = c(1,1)
            )
            )
            
            tmt_intensity_cols = get_tmt_intensity_columns(dt, des)
            expect_equal(tmt_intensity_cols,
                         c("reporter intensity corrected 1  EV",
                           "reporter intensity corrected 2  EV",
                           "reporter intensity corrected 3  EV",
                           "reporter intensity corrected 4  EV", 
                           "reporter intensity corrected 5  EV"))
          })

test_that("Returns tmt intensity when not corrected is there and no file name", 
          {
            
            des = as.data.table(list(
              condition = c("A","B","A", "B","A"),
              fraction = c(0,0,0,0,0,0,0),
              experiment = c("","","","",""),
              remove = c(FALSE,FALSE,FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,1,2,2,3),
              reporter_channel = c(1,2,3,4,5)
            ))
            
            dt = as.data.table(list(
              fasta_headers = c("blablabla","blablabla2"),
              `reporter intensity not corrected 1` = c(1,1),
              `reporter intensity not corrected 2` = c(1,1),
              `reporter intensity not corrected 3` = c(1,1),
              `reporter intensity not corrected 4` = c(1,1),
              `reporter intensity not corrected 5` = c(1,1)
            )
            )
            
            tmt_intensity_cols = get_tmt_intensity_columns(dt, des)
            expect_equal(tmt_intensity_cols,
                         c("reporter intensity not corrected 1",
                           "reporter intensity not corrected 2",
                           "reporter intensity not corrected 3",
                           "reporter intensity not corrected 4", 
                           "reporter intensity not corrected 5"))
          })

test_that("Returns tmt intensity when no file name", 
          {
            
            des = as.data.table(list(
              condition = c("A","B","A", "B","A"),
              fraction = c(0,0,0,0,0,0,0),
              experiment = c("","","","",""),
              remove = c(FALSE,FALSE,FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,1,2,2,3),
              reporter_channel = c(1,2,3,4,5)
            ))
            
            dt = as.data.table(list(
              fasta_headers = c("blablabla","blablabla2"),
              `reporter intensity corrected 1` = c(1,1),
              `reporter intensity corrected 2` = c(1,1),
              `reporter intensity corrected 3` = c(1,1),
              `reporter intensity corrected 4` = c(1,1),
              `reporter intensity corrected 5` = c(1,1)
            )
            )
            
            tmt_intensity_cols = get_tmt_intensity_columns(dt, des)
            expect_equal(tmt_intensity_cols,
                         c("reporter intensity corrected 1",
                           "reporter intensity corrected 2",
                           "reporter intensity corrected 3",
                           "reporter intensity corrected 4", 
                           "reporter intensity corrected 5"))
          })


