# Similar to reformat_rna.R but for sections marked ***
rm(list = ls())
library(tidyverse)
library(Biobase) 

########################################################
# LOAD DATA
########################################################
########################################################
# Same data files, different paths/servers
########################################################
exec_from_uzh  <-  TRUE
if(exec_from_uzh){
  
  
  # phenotype from lauren, covariates: see email.
  ph = "/Volumes/Data/Addhealth/RNA/180626_RNAControls.dta" %>%
    haven::read_dta()   
  ph = ph %>% mutate(AID = aid) 
  # preprocessed expression matrix
  ex_preproc     = "/Volumes/Data/Addhealth/RNA/AddHealthYr1 Plates1-12 - Gene CPM Log2.txt" %>% 
    read.table(stringsAsFactors = FALSE, header = TRUE) 
  # raw expression
  ex  = "/Volumes/Data/Addhealth/RNA/AddHealthYr1 Plates1-12 - ReadsPerGene ENSG Raw.txt" %>% 
    read.table(stringsAsFactors = FALSE, header = TRUE) 
  # relate expression sample ids - addhealth aids
  ids = "/Volumes/Data/Addhealth/RNA/AID_VialID.csv" %>%
    read.csv(stringsAsFactors = FALSE, header = TRUE) #%>% 
  ids = ids %>% mutate(VialID = str_c("X",VialID ), AID = as.character(AID))
  # RNA quality metrics
  quality <- read.table("/volumes/data/addhealth/RNA/AddHealthYr1_Plates1-12_Log2CPM_QCMetrics.txt",  header = TRUE, dec = ".") %>% as_tibble
  quality = quality %>% 
    mutate(VialID = str_c("X",Sample), 
           Plate = factor(Plate), 
           AvgCorrelogram100 = as.numeric(AvgCorrelogram100)) %>% 
    select(VialID, Plate, AvgCorrelogram100)
  
} else {
  
  # phenotype from lauren, covariates: see email.
  ph = "/ifs/sec/cpc/addhealth/users/lgaydosh/180626_RNAControls.dta" %>%
    haven::read_dta()   
  ph = ph %>% mutate(AID = aid) 
  # preprocessed expression matrix
  ex_preproc     = "/ifs/sec/cpc/addhealth/RNA/AddHealthYr1 Plates1-12 - Gene CPM Log2.txt" %>% 
    read.table(stringsAsFactors = FALSE, header = TRUE) 
  # raw expression
  ex  = "/ifs/sec/cpc/addhealth/RNA/AddHealthYr1 Plates1-12 - ReadsPerGene ENSG Raw.txt" %>% 
    read.table(stringsAsFactors = FALSE, header = TRUE) 
  # relate expression sample ids - addhealth aids
  ids = "/ifs/sec/cpc/addhealth/RNA/IDs/AID_VialID.csv" %>%
    read.csv(stringsAsFactors = FALSE, header = TRUE) #%>% 
  ids = ids %>% mutate(VialID = str_c("X",VialID ), AID = as.character(AID))
  
}

# *** assumption checked below ***
# write file in order to query http://www.ensembl.org
ensgcode= ex$Gene # for posterity
ex$Gene = ex_preproc$Gene
write.csv(ensgcode, "data/ensgcode.csv") 

# check above assumption mapping ensembl to hgnc
if(0){
  mapping_used_above = tibble(ensembl_gene_id = ex$Gene, hgnc_symbol = ex_preproc$Gene)  
  
  # find mapping on biomart, for confirmation
  library(biomaRt)
  my_mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
  my_arg_var      = c("ensembl_gene_id")       # argument variable of my search function
  my_arg_var_val  = list(ex$Gene)              # argument value of my search function
  my_return_vars  = c("hgnc_symbol", "entrezgene", my_arg_var) # return value of my search function
  (hits=getBM(attributes=my_return_vars,
              filters=my_arg_var,
              values=my_arg_var_val,
              mart=my_mart))
  
  # check the last two collumns are equal
  hits %>% left_join(mapping_used_above, by = "ensembl_gene_id") %>% head 
}

# remove duplicate genes
# fix ex rownames
ex = ex[!duplicated(ex$Gene), ]
rownames(ex) = ex$Gene[!duplicated(ex$Gene)]
# fix colnames
ex = ex %>% select(-Gene) #  your job is done

# primary key = VialID (commented out, I could join in Laurens phenotype data)
ph = ids %>% left_join(ph, by = "AID")

# select those augmented phenotype matrix with corresponding samples
sampleIDs = colnames(ex) 
rownames(ph) = ph$VialID # needed for indexing (next line)
ph = ph[sampleIDs, ] %>% left_join(quality, by = "VialID") 
rownames(ph) = ph$VialID # needed by ExpressionSet() function

# make an expression set
# check 
all.equal(colnames(ex), rownames(ph)) 
phenData = new("AnnotatedDataFrame", data = ph) 
dat = ExpressionSet(assayData = ex %>% as.matrix, 
                    phenoData = phenData)  

# Identify  11 duplicated subjects seperately 
dupes = 
  pData(dat)$AID %>% 
  fct_count() %>% 
  filter(n == 2) %>% 
  mutate(f = as.character(f)) %>% 
  select(f) %>% 
  unlist 

dat_dupes = dat[, dat$AID %in% dupes]
# save ***
dat_dupes %>% saveRDS(file = "data/dt_dupes.rds")

# final dataset: remove one of the above duplicates (the second)
dat = dat[, !duplicated(dat$AID) ] 
# save  ***
dat %>% saveRDS(file = "data/dt.rds")

# for convenience
eD = dat %>% exprs 
pD = dat %>% pData 
fD = dat %>% featureData 
aD = cbind(pD,t(eD)) %>% as_tibble 
# all data ***
aD %>% saveRDS(file = "data/dtTidy.rds")
