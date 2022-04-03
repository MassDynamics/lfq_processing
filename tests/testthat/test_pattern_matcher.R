library(LFQProcessing)
library(testthat)
library(tools)


test_that("Match Uniprot correctly", 
          {
          
            headerString = 'tr|Q8N2H2|Q8N2H2_HUMAN cDNA FLJ90785 fis, clone THYRO1001457, moderately similar to H.sapiens protein kinase C mu OS=Homo sapiens OX=9606 PE=2 SV=1'
            expect_true(
              match_pattern_name(headerString) == "uniprot"
            )  
          }
)



test_that("Match ups accession correctly", 
          {
            
            headerString = 'P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human);'
            expect_true(
              match_pattern_name(headerString) == "accession"
            )  
          }
)

test_that("Match Accession correctly", 
          {
            
            headerString = 'tr|Q8N2H2|Q8N2H2_HUMAN'
            expect_true(
              match_pattern_name(headerString) == "accession"
            )  
          }
)

test_that("Match con correctly", 
          {
            
            headerString = 'CON__ENSEMBL:ENSBTAP00000038253'
            expect_equal(
              match_pattern_name(headerString), "con"
            )  
          }
)

test_that("Match hprr correctly", 
          {
            
            headerString = 'HPRR670317'
            expect_true(
              match_pattern_name(headerString) == "hppr"
            )  
          }
)


test_that("Match missing correctly", 
          {
            
            headerString = ''
            expect_true(
              match_pattern_name(headerString) == "missing"
            )  
          }
)


test_that("Match null as missing", 
          {
            
            headerString = NULL
            expect_true(
              match_pattern_name(headerString) == "missing"
            )  
          }
)
          
test_that("Match na as missing", 
          {
            
            headerString = NA
            expect_true(
              match_pattern_name(headerString) == "missing"
            )  
          }
)

test_that("Match semicolons as missing", 
          {
            
            headerString = ";;;"
            expect_true(
              match_pattern_name(headerString) == "missing"
            )  
          }
)

test_that("Match semicolon as missing", 
          {
            
            headerString = ";"
            expect_true(
              match_pattern_name(headerString) == "missing"
            )  
          }
)