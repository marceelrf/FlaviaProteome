#https://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html


library(tidyverse)
library(readxl)
library(limma)
library(qvalue)

dat <- read.csv("http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/TMT_experiment.csv") 
raw_data <- read_xlsx(path = "Data/Protein_TMT_Quant.xlsx")


source("http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/read.peptides.r")
source("http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/quantify.proteins.r")
source("http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/eb.fit.r")


(cha <- colnames(raw_data[17:32]))

dat <- read.peptides(raw_data, cha) 
dim(dat) # 19047 peptides


# Create a BIOBASE obj ----------------------------------------------------

exp_mat <- raw_data[,17:32]

rownames(exp_mat) <- raw_data$PG.ProteinAccessions


feature_data <- raw_data[,1:16]
rownames(feature_data) <- raw_data$PG.ProteinAccessions

pheno_data <- read_xlsx(path = "Data/Sample-Info.xlsx")
rownames(pheno_data) <- pheno_data$BioReplicate

pheno_data <- pheno_data[order(rownames(pheno_data)),]
rownames(pheno_data) <- pheno_data$BioReplicate


expset <-
  Biobase::ExpressionSet(assayData = as.matrix(exp_mat),
                       phenoData = Biobase::AnnotatedDataFrame(pheno_data),
                       featureData = Biobase::AnnotatedDataFrame(feature_data))
