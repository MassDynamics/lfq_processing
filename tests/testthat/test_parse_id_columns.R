library(LFQProcessing)
library(testthat)
library(tools)


test_that("Heterogenous cases runs and produces correct output", 
          {
            prot <- as.data.table(list(id = c(1,2,3,4,5, 6),
                                       `majority protein ids` = c("P0DPI2;A0A0B4J2D5", "A0A0U1RRL7", "A0AV96----", "HPRR670317", "", "junk"),
                                       `fasta headers` = c(
                                         "sp|P0DPI2|GAL3A_HUMAN Glutamine amidotransferase-like class 1 domain-containing protein 3A, mitochondrial OS=Homo sapiens OX=9606 GN=GATD3A PE=1 SV=1;sp|A0A0B4J2D5|GAL3B_HUMAN Glutamine amidotransferase-like class 1 domain-containing protein 3B, mitochondr", 
                                         "sp|A0A0U1RRL7|MMPOS_HUMAN", 
                                         ";",
                                         "HPRR670317",
                                         "P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human);", "junk")))
            
          prot <- parse_id_columns(prot)  
          
          expect_true("Accession_id" %in% colnames(prot))
          expect_true("Gene" %in% colnames(prot))
          expect_true("ProteinDescription" %in% colnames(prot))
          
          expect_equal(prot$Gene[[1]], "GATD3A")
          expect_equal(prot$Accession_id, c("P0DPI2", "A0A0U1RRL7", "A0AV96----", "HPRR670317", "P02768",     "junk" ))
          expect_equal(prot$ProteinDescription, c(" Glutamine amidotransferase-like class 1 domain-containing protein 3A, mitochondrial", 
                                                  "sp|A0A0U1RRL7|MMPOS_HUMAN", "", "HPRR670317", "P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human)",     "junk" ))
          }
          
          )

test_that("Uniprot cases handled well", 
          {
            prot <- as.data.table(list(id = c(1,2,3,4,5,6),
                                       `majority protein ids` = c("","","","","",""),
                                       `fasta headers` = c(
                                         'sp|Q8I6R7|ACN2_ACAGO Acanthoscurrin-2 (Fragment) OS=Acanthoscurria gomesiana OX=115339 GN=acantho2 PE=1 SV=1',
                                         'sp|P27748|ACOX_CUPNH Acetoin catabolism protein X OS=Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337) OX=381666 GN=acoX PE=4 SV=2',
                                         'sp|P04224|HA22_MOUSE H-2 class II histocompatibility antigen, E-K alpha chain OS=Mus musculus OX=10090 PE=1 SV=1',
                                         'tr|Q3SA23|Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1',
                                         'tr|Q8N2H2|Q8N2H2_HUMAN cDNA FLJ90785 fis, clone THYRO1001457, moderately similar to H.sapiens protein kinase C mu OS=Homo sapiens OX=9606 PE=2 SV=1',
                                         '>sp|A0A075B6I9|LV746_HUMAN Immunoglobulin lambda variable 7-46 OS=Homo sapiens OX=9606 GN=IGLV7-46 PE=3 SV=4')
            ))
            prot <- parse_id_columns(prot)

            expect_equal(unique(unlist(lapply(prot$`fasta headers`, match_pattern_name))), "uniprot")            
            expect_equal(prot$Accession_id,c("Q8I6R7", "P27748","P04224","Q3SA23","Q8N2H2","A0A075B6I9"))
            expect_equal(prot$Gene[[1]],"acantho2")
            expect_equal(prot$Gene[[2]],"acoX")
            expect_equal(prot$Gene[[4]],"nef")
            expect_equal(prot$Gene[[6]],"IGLV7-46")
          }
          
)

test_that("Uniprot Isoform cases handled well", 
          {
            prot <- as.data.table(list(id = c(1,2),
                                       `majority protein ids` = c("",""),
                                       `fasta headers` = c(
                                         'sp|P52480-2|KPYM_MOUSE Isoform M1 of Pyruvate kinase PKM OS=Mus musculus GN=Pkm',
                                         'sp|P52480|KPYM_MOUSE Pyruvate kinase PKM OS=Mus musculus GN=Pkm PE=1 SV=4')
            ))
            prot <- parse_id_columns(prot)
            
            expect_equal(unique(unlist(lapply(prot$`fasta headers`, match_pattern_name))), "uniprot")            
            expect_equal(prot$Accession_id,c("P52480-2", "P52480"))
            expect_equal(prot$Gene[[1]],"Pkm")
            expect_equal(prot$Gene[[2]],"Pkm")
            expect_equal(prot$ProteinDescription[[1]]," Isoform M1 of Pyruvate kinase PKM")
            expect_equal(prot$ProteinDescription[[2]]," Pyruvate kinase PKM")
          }
          
)

test_that("HPPR cases handled well", 
          {
            prot <- as.data.table(list(id = c(1,2,3,4,5,6,7),
                                       `majority protein ids` = c("","","","","","",""),
                                       `fasta headers` = c(
                                         "HPRR1010001",
                                         "HPRR1310038",
                                         "HPRR1320009",
                                         "HPRR1370013",
                                         "HPRR1370116",
                                         "HPRR140662",
                                         "HPRR140746")
            ))
            prot <- parse_id_columns(prot)
            
            expect_equal(unique(unlist(lapply(prot$`fasta headers`, match_pattern_name))), "hppr")     
            expect_equal(prot$Accession_id, c(
              "HPRR1010001",
              "HPRR1310038",
              "HPRR1320009",
              "HPRR1370013",
              "HPRR1370116",
              "HPRR140662",
              "HPRR140746"))       
          }
          
)

test_that("Missing cases handled well",
          {
            prot <- as.data.table(list(id = c(1,2,3,4),
                                       `majority protein ids` = c("P08397",
                                                                  "P08559;P29803",
                                                                  "P08603",
                                                                  "P08670;P41219;P17661"),
                                       `fasta headers` = c("",
                                                           ";",
                                                           "",
                                                           ";;")
            ))
            prot <- parse_id_columns(prot)
            
            expect_equal(unique(unlist(lapply(prot$`fasta headers`, match_pattern_name))), "missing")     
            expect_equal(prot$Accession_id, c("P08397", "P08559", "P08603", "P08670"))
          })


test_that("Test modified accession",
          {
            prot <- as.data.table(list(id = c(1,2),
                                       `majority protein ids` = c("QPatient53;QPatient52;Q14498", "QPatient53;QPatient52"),
                                       `fasta headers` = c("sp|QPatient53|RBM39_HUMAN Isoform 3 of RNA-binding protein 39 OS=Homo sapiens OX=9606 GN=RBM39;sp|QPatient52|RBM39_HUMAN Isoform 2 of RNA-binding protein 39 OS=Homo sapiens OX=9606 GN=RBM39;sp|Q14498|RBM39_HUMAN RNA-binding protein 39 OS=Homo sapiens OX=9606 G",
                                                           "sp|QPatient53|RBM39_HUMAN Isoform 3 of RNA-binding protein 39 OS=Homo sapiens OX=9606 GN=RBM39;sp|QPatient52|RBM39_HUMAN Isoform 2 of RNA-binding protein 39 OS=Homo sapiens OX=9606 GN=RBM39")
            ))
            prot <- parse_id_columns(prot)
            
            expect_equal(unique(unlist(lapply(prot$`fasta headers`, match_pattern_name))), c("uniprot"))    
            expect_equal(prot$Accession_id, c("QPatient53","QPatient53"))
          })


test_that("Test null id",
          {
            prot <- as.data.table(list(id = c(1),
                                       `majority protein ids` = c(""),
                                       `fasta headers` = c("Not Null")
            ))
            prot <- parse_id_columns(prot)

            expect_equal(prot$Accession_id, c("Not Null"))
          })

