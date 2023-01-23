## Kinetic analysis of equine mesenchymal stromal cell secreted proteins
## Stable Isotope Dynamic Labelling of Secretomes (SIDLS)

# load packages
library(stringr)
library(dplyr)
library(tidyr)

# load MaxQuant 'evidence' file with peptide-level data
evi <- read.csv('C:/Users/aggiej/Documents/Proteomics/SILAC_0421/Evidence_no_requant.csv')

# remove contaminants and reverse hits flagged by MaxQuant
evi0 <- evi %>% filter(Potential.contaminant != '+') %>% filter(Reverse != '+')

## Add columns with time point and RIA for each peptide in evidence file
evi1 <- evi0 %>% mutate(time = str_remove(Experiment, '[H]'), RIA = Intensity.H/(Intensity.H + Intensity.L))

# filter peptides identified in at least 3/4 time points
evi2 <- evi1 %>% group_by(Sequence) %>% summarise(count = n_distinct(Experiment)) %>% filter(count >=3)

# identified proteins with peptide data identified in at least 3/4 time points
evi3 <- evi2 %>% left_join(evi1, by = 'Sequence')
list_proteins <- unique(evi3$Leading.razor.protein)

# fitting optimal curve to represent RIA change over time (labelling trajectory)
# define function
fn_regresults <- function(iprot) {
  
  tryCatch({
    
    # For protein select RIA + time
    evi1_subset <- evi1 %>% 
      filter(Leading.proteins == iprot) %>% 
      mutate(time = as.numeric(time, levels = c(1, 2, 6, 24))) %>%
      select(RIA,time) %>%
      na.omit()
    
    ## Fit model
    nls1 <- nls(RIA ~ (1 - exp(-k*time)),
                data = evi1_subset,
                start = list(k = 0.2))
    
    return(
      data.frame(
        prot = iprot,
        broom::tidy(nls1)
      )
    )
    
  }, error=function(e){})
}

    
    

# apply function to proteins with peptides identified at 3/4 time points
# estimate = k value
out1 <- purrr::map_dfr(list_proteins, ~fn_regresults(.))

# adjust p-values for multiple comparisons
out2 <- out1 %>% mutate(q.value = p.adjust(out1$p.value, method = 'BH'))

# calculating difference in total protein abundance between 24 and 6h time point (P)
# define function
fn_abundprot <- function (iprot) {
  
  # extract summed peak intensity for each peptide
  abund <- evi1 %>% 
    filter(Leading.proteins == iprot) %>% 
    mutate(time = as.numeric(time, levels = c(1, 2, 6, 24))) %>%
    select(Intensity, time) %>%
    na.omit()
  
  # calculate total mean intensity for each time point
  mean_abund <- abund %>% 
    group_by(time) %>%
    summarise(mean = mean(Intensity), SD = sd(Intensity))
  
  # select relevant columns and change data format from long to wide
  mean_abund1 <- mean_abund %>% 
    select (time, mean) %>% 
    pivot_wider( names_from = time, values_from = mean)
  
  # calculate abundance over last 18h
  P = mean_abund1$`24` - mean_abund1$`6`
  P
  
}

# apply function to list of proteins previously characterised with k value
list <- unique(out1$prot)
out3 <- lapply(list, fn_abundprot)
out4 <- data.frame(Reduce(rbind, out3))
colnames(out4) <- 'P'

# join RIA trajectory data and P value
out5 <- out1 %>% cbind(out4)

# calculate flux of protein from intra- to extracellular space
out6 <- out5 %>% mutate(flux = estimate * P)

# calculate log102flux

out7 <- out6 %>% mutate(log10flux = log10(flux))

# save results

write.csv(out7, 'Supplementary_File_2.csv')
