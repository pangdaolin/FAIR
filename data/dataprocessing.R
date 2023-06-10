library(plyr)
library(caret)
library(openxlsx)
library(readr)

lc_study_otu_table <- read_tsv("/home/daolinpang/project1/data/lc_study_otu_table.tsv")
lc_study_mapping_file <- read_tsv("/home/daolinpang/project1/data/lc_study_mapping_file.tsv")
data.id <- read.xlsx("/home/daolinpang/project1/data/pbio.2003862.s004.xlsx")	
imp.otu <- read.xlsx("/home/daolinpang/project1/data/pbio.2003862.s040.xlsx")
