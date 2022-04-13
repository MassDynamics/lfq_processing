library(LFQProcessing)
library(testthat)
library(tools)

test_that("Score exists and doesn't replace standard fasta headers", 
          {
            
            prot <- as.data.table(list(id = c(1,2,3,4,5, 6),
                                       `majority protein ids` = c("P0DPI2;A0A0B4J2D5", "A0A0U1RRL7", "A0AV96----", "HPRR670317", "", "junk"),
                                       `fasta headers` = c(
                                         "sp|P0DPI2|GAL3A_HUMAN Glutamine amidotransferase-like class 1 domain-containing protein 3A, mitochondrial OS=Homo sapiens OX=9606 GN=GATD3A PE=1 SV=1;sp|A0A0B4J2D5|GAL3B_HUMAN Glutamine amidotransferase-like class 1 domain-containing protein 3B, mitochondr", 
                                         "sp|A0A0U1RRL7|MMPOS_HUMAN", 
                                         ";",
                                         "HPRR670317", "A", "B")))
            
            for (col in c("protein ids", "q-value","peptide counts (all)", "peptide counts (razor+unique)", 
                          "peptide counts (unique)", "number of proteins", "peptides", "razor + unique peptides", 
                          "unique peptides", "sequence coverage [%]", "unique + razor sequence coverage [%]", 
                          "unique sequence coverage [%]", "mol. weight [kda]" , 'sequence length', "peptide ids", 
                          'evidence ids', 'ms/ms ids', 'best ms/ms', 
                          'lfq intensity sample1rep1', 'lfq intensity sample1rep2','lfq intensity sample1rep3')){
              prot[,col] = 1
            }
            
            des = as.data.table(list(
              file_name = c("A","B","C"),
              fraction = c(0,0,0),
              experiment = c("bdcya","bdcya","bdcya"),
              remove = c(FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,2,3)
            ))
            
            prot <- lfq_proteinGroup_txt_formatter(prot, des)
            
            expect_true("score" %in% colnames(prot))
            expect_true("lfq intensity sample1rep1" %in% colnames(prot))
            expect_true("lfq intensity sample1rep2" %in% colnames(prot))
            expect_true("lfq intensity sample1rep3" %in% colnames(prot))
            
            expect_equal(prot$`fasta headers`, c(
              "sp|P0DPI2|GAL3A_HUMAN Glutamine amidotransferase-like class 1 domain-containing protein 3A, mitochondrial OS=Homo sapiens OX=9606 GN=GATD3A PE=1 SV=1;sp|A0A0B4J2D5|GAL3B_HUMAN Glutamine amidotransferase-like class 1 domain-containing protein 3B, mitochondr", 
              "sp|A0A0U1RRL7|MMPOS_HUMAN", 
              ";",
              "HPRR670317", "A","B"))
          }
)

test_that("Replace empty fasta headers test", 
          {
            
            prot <- as.data.table(list(id = c(1,2,3,4,5, 6),
                                       `majority protein ids` = c("P0DPI2;A0A0B4J2D5", "A0A0U1RRL7", "A0AV96----", "HPRR670317", "", "junk"),
                                       `fasta headers` = c(
                                         "", 
                                         "", 
                                         ";",
                                         "",
                                         "P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human);", "junk")))
            
            for (col in c("protein ids", "q-value","peptide counts (all)", "peptide counts (razor+unique)", 
                          "peptide counts (unique)", "number of proteins", "peptides", "razor + unique peptides", 
                          "unique peptides", "sequence coverage [%]", "unique + razor sequence coverage [%]", 
                          "unique sequence coverage [%]", "mol. weight [kda]" , 'sequence length', "peptide ids", 
                          'evidence ids', 'ms/ms ids', 'best ms/ms', 
                          'lfq intensity sample1rep1', 'lfq intensity sample1rep2','lfq intensity sample1rep3')){
              prot[,col] = 1
            }
            
            des = as.data.table(list(
              file_name = c("A","B","C"),
              fraction = c(0,0,0),
              experiment = c("bdcya","bdcya","bdcya"),
              remove = c(FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,2,3)
            ))
            
            prot <- lfq_proteinGroup_txt_formatter(prot, des)
            
            expect_true("score" %in% colnames(prot))
            expect_true("lfq intensity sample1rep1" %in% colnames(prot))
            expect_true("lfq intensity sample1rep2" %in% colnames(prot))
            expect_true("lfq intensity sample1rep3" %in% colnames(prot))
            
            expect_equal(sum(duplicated(prot$`fasta headers`)),0)
          }
)

test_that("Remove empty id's", 
          {
            
            prot <- as.data.table(list(id = c(1,2,3,4, NA, "", NULL),
                                       `majority protein ids` = c("P0DPI2;A0A0B4J2D5", "A0A0U1RRL7", "A0AV96----", "HPRR670317", "", "junk"),
                                       `fasta headers` = c(
                                         "", 
                                         "", 
                                         ";",
                                         "",
                                         "P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human);", "junk")))
            
            for (col in c("protein ids", "q-value","peptide counts (all)", "peptide counts (razor+unique)", 
                          "peptide counts (unique)", "number of proteins", "peptides", "razor + unique peptides", 
                          "unique peptides", "sequence coverage [%]", "unique + razor sequence coverage [%]", 
                          "unique sequence coverage [%]", "mol. weight [kda]" , 'sequence length', "peptide ids", 
                          'evidence ids', 'ms/ms ids', 'best ms/ms', 
                          'lfq intensity sample1rep1', 'lfq intensity sample1rep2','lfq intensity sample1rep3')){
              prot[,col] = 1
            }
            
            des = as.data.table(list(
              file_name = c("A","B","C"),
              fraction = c(0,0,0),
              experiment = c("bdcya","bdcya","bdcya"),
              remove = c(FALSE,FALSE,FALSE),
              mqExperiment = c("sample1rep1","sample1rep2","sample1rep3"),
              Replicate = c(1,2,3)
            ))
            
            prot <- lfq_proteinGroup_txt_formatter(prot, des)
            
            expect_true("score" %in% colnames(prot))
            expect_true("lfq intensity sample1rep1" %in% colnames(prot))
            expect_true("lfq intensity sample1rep2" %in% colnames(prot))
            expect_true("lfq intensity sample1rep3" %in% colnames(prot))
            
            expect_equal(nrow(prot),4)
          }
)