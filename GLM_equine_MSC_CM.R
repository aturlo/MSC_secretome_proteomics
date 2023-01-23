## Analysis of secreted proteins in equine mesenchymal stromal cell conditioned media using generalised linear models (GLMs)
## dataset obtained through label-free mass spectrometry proteomics

# load packages
library(broom)
library('tidyr')
library("tidyverse")
library("ggplot2")
library("limma")

# load label-free quantification (LFQ) protein-level data 

all <- read.csv('C:/Users/aggiej/Documents/Proteomics/MSC_secretome_study/MaxQuant_output/All_LFQ_filtered.csv')

# load subset of proteins identified as secreted with SIDLS
sec <- read.csv('C:/Users/aggiej/Documents/Proteomics/SILAC_0421/Kinetic analysis/Subsets/SILAC_BMSC_secretome.csv')

# format and merge datasets to obtain LFQ values only for the proteins of interest
sec1 <- as.vector(sec$Protein)
all_long <- all %>% pivot_longer(!c('X', 'Age', 'Tissue', 'Sex'), 
                                 names_to = 'Protein', 
                                 values_to = 'LFQintensity')

datafin2 <- all_long %>% filter(Protein %in% sec1)
unique(datafin2$Protein)

# Apply GLM with Gamma distribution and log link to test effect of donor age, sex and tissue source on LFQ intensity
# define function

reg_coefs <- function(prot){
  subset <- datafin2 %>%
    filter(Protein == prot) %>%
    drop_na()
  
  model <- glm(LFQintensity ~ Age + Sex + Tissue, 
               family = Gamma(link = 'log'), 
               data = subset)
  
  coefs <- as.data.frame(coef(summary(model)))
  coefs <- coefs %>% mutate(predictor = rownames(coefs))
  
}

# extract list of target protein names
prots <- unique(datafin2$Protein)

# apply function to the list of target proteins
coefs_gamma <- lapply(prots, reg_coefs)

# transform the list of data frames to a data frame
coefs_gamma1 <- bind_rows(coefs_gamma, .id = 'prt') %>%
  mutate(protein = rep(prots, each = 4)) %>%
  remove_rownames() %>%
  rename(p.value = 'Pr(>|t|)')

# apply Benjamini-Hochberg correction for multiple comparisons

coefs_gamma2 <- coefs_gamma1 %>% filter(predictor != '(Intercept)') # exclude p-values associated with intercept
coefs_gamma3 <- coefs_gamma2 %>% mutate(adjusted_p.value = p.adjust(p.value, method = 'BH'))


##  Log2 fold change calculation
# log-transform LFQ intensity values 
datafin2_log <- datafin2 %>% mutate(LFQintensity_log = log2(LFQintensity))


# calculate log fold change between samples from different tissue types
fc_bm <- datafin2_log %>%
  group_by(Protein) %>%
  filter(Tissue == 'Bone marrow') %>%
  mutate(mean_BM = mean(LFQintensity_log, na.rm = TRUE)) %>%
  filter(!duplicated(Protein)) %>%
  select(Protein, mean_BM)


fc_adi <- datafin2_log %>%
  group_by(Protein) %>%
  filter(Tissue == 'Adipose') %>%
  mutate(mean_adi = mean(LFQintensity_log, na.rm = TRUE)) %>%
  filter(!duplicated(Protein)) %>%
  select(Protein, mean_adi)


fc <- full_join(fc_bm, fc_adi, by = 'Protein')
fc <- fc %>% mutate(log_fc = mean_BM - mean_adi)


# calculate log fold change between samples from donors of different sexes

fc_m <- datafin2_log %>%
  group_by(Protein) %>%
  filter(Sex == 'M') %>%
  mutate(mean_M = mean(LFQintensity_log, na.rm = TRUE)) %>%
  filter(!duplicated(Protein)) %>%
  select(Protein, mean_M)

fc_f <- datafin2_log %>%
  group_by(Protein) %>%
  filter(Sex == 'F') %>%
  mutate(mean_F = mean(LFQintensity_log, na.rm = TRUE)) %>%
  filter(!duplicated(Protein)) %>%
  select(Protein, mean_F)

fcs <- full_join(fc_m, fc_f, by = 'Protein')
fcs <- fcs %>% mutate(log_fc = mean_M - mean_F)


# join results of GLM and log2 fold change calculation 
fc_glm <- coefs_gamma3 %>% full_join(fc, by = c('protein' = 'Protein')) %>% full_join(fcs, by = c('protein' = 'Protein') )
fc_glm1 <- fc_glm %>% mutate(log_FC = ifelse(predictor == 'TissueBone marrow', log_fc.x, 
                                             ifelse(predictor == 'SexM', log_fc.y, NA)))

# change column names
fc_glm2 <- fc_glm1 %>% rename(GLM_coefficient = Estimate,
                              GLM_std.error = 'Std. Error',
                              GLM_t.value = 't value',
                              GLM_p.value = p.value,
                              GLM_adjusted_p.value = adjusted_p.value,
                              Protein_accession_no = protein)

fc_glm3 <- fc_glm2 %>% select(c(GLM_coefficient, GLM_std.error, GLM_t.value, GLM_p.value, GLM_adjusted_p.value, predictor, Protein_accession_no, log_FC))

# change column order
fc_glm3 <- fc_glm3[, c(7, 6, 1, 2, 3, 4, 5, 8)]

# save results
write.csv(fc_glm3, 'Supplementary_File_4.csv')
