## Processing of label-free mass spectrometry data

# load packages

library('tidyr')
library("tidyverse")
library("ggplot2")

# load data (MaxQuant 'proteinGroups' file)
data <- read.csv('C:/Users/aggiej/Documents/Proteomics/MSC_secretome_study/MaxQuant_output/proteinGroups.csv')

# load metadata
meta <- read.csv('C:/Users/aggiej/Documents/Proteomics/MSC_secretome_study/MaxQuant_output/MSC_metadata.csv')

# filter based on 2 unique peptides/protein group

data2up <- data %>% filter(Unique.peptides >=2)

# subset variables to keep only useful ones
data2up <- data2up %>% select("Majority.protein.IDs",  "Fasta.headers", "Q.value", "Score", contains("LFQ.intensity"),  # create vector with variable names
                              "Reverse", "Potential.contaminant")

# clean up quantitative column names to leave only sample id
colnames(data2up) <- sub('LFQ.intensity.',
                         '',
                         colnames(data2up))

## MISSING VALUES
# replace 0 in MaxQuant output with NA
data2up[data2up == 0] <- NA
data2up$Q.value[is.na(data2up$Q.value)] <- 0 # apart from the q value column

# filter contaminants and reverse hits
dataClean <- data2up %>% filter(Reverse == '') %>% filter(Potential.contaminant == '')

# filter q value < 0.01
dataClean <- dataClean %>% filter(Q.value <= 0.01)

## FORMAT DATA FOR STATISTICAL ANALYSIS
# subset dataframe with protein ID and quantitative columns only
stat <- dataClean %>% select(1, 5:36)

# split majority protein names into columns
statX <- str_split_fixed(stat$Majority.protein.IDs, ';', 4)
statX <- statX[,1] # keep only 1 first protein
statX <- as.data.frame(statX) # convert vector into data frame
statY <- cbind(stat, statX) # add column with first majority protein id to quantitative data
stat2 <- statY %>% select(statX, everything()) # reorder last column as the 1st
stat2 <- stat2[,-2] # remove the column with multiple protein ids


# transpose data frame
stat3 <- t(stat2)
stat3 <- as.data.frame(stat3) # save matrix as data frame
colnames(stat3) <- stat3[1,] # assign first row as column names
stat3 <- stat3[-1,] # remove first row

# add metadata
predictors <- meta[,c(1,3,4)]
stat4 <- cbind(stat3, predictors)

# remove samples with missing meta data
stat5 <- stat4[-c(17, 21),]

# save to csv file
write.csv(stat5, 'All_LFQ_filtered.csv')
