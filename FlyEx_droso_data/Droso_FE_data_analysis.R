# -----------------------------------------
# Droso_FE_data_analysis. R
#   
#   perform clean up and EDA for FlyEx Drosophila gap gene expression data
#   then output cleaned data for other uses
#   database: FlyEx
# -----------------------------------------

# -----------------------------------------
# Author:Jianhong Chen
# Date: 04-23-2019
# -----------------------------------------


# -----------------------------------------------------------------------------------------------------
# load require library
library('ggplot2')
library('tidyverse') # inlcude "tidyr", "readr", "dplyr", "purrr"
library('reshape2') 

# -----------------------------------------------------------------------------------------------------
# intial set up
ls() # list all objects 
gc() # force R to release memory of deleted variables 
rm(list = ls()) # delete all variables in the environment
set.seed(78) # set random seed
# -----------------------------------------------------------------------------------------------------

# converting into standarized AP position (100 points)
Strip_Data = function(x, sep,  gen1, gen2){
  # bin the intensity measurements of the gap genes into specified strips by taking the mean value:
  #   inputs: x = data file;
  #           sep = width of the strip
  #           gen1, gen2 = interested genes
  #   output: a dataframe of the binned (averaged) of the gap gene 
  oneAvg = data.frame()
  for (n in seq(0, (100-sep),sep)){
    one = subset(x, (x$AP > n) & (x$AP < (n+sep)))
    oneAvg = rbind(oneAvg, c(mean(one[[gen1]]), mean(one[[gen2]])))
  }
  colnames(oneAvg) = c(gen1, gen2)
  return(oneAvg)
}


# -----------------------------------------------------------------------------------------------------
#  1. load the FlyEx data files
# -----------------------------------------------------------------------------------------------------
# data for gt & kni:
dir1 = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Gap_Gene_of_Drosophila/FlyEx_data/gt_kni14At1-8_w_R3/txt/byEmbryos"
setwd(dir1)
f1 = list.files()
gt_kni_raw = map(f1, read.table, header = T, col.names = c('Nucleus_Num', 'AP', 'DV', 'eve', 'kni', 'gt'))

# data for kr & hb:
dir2 = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Gap_Gene_of_Drosophila/FlyEx_data/kr_hb14At1-8_w_R3/txt/byEmbryos"
setwd(dir2)
f2 = list.files()
kr_hb_raw = map(f2, read.table, header = T, col.names = c('Nucleus_Num', 'AP', 'DV', 'eve', 'kr', 'hb'))

# quick plot of the raw data without any processing
gt_kni_exp = data.frame()
kr_hb_exp = data.frame()
for (i in 1:8) {
  gt_kni_temp = gt_kni_raw[[i]] %>%
    mutate('time' = rep(i, dim(gt_kni_raw[[i]])[1]))
  gt_kni_exp = rbind(gt_kni_exp, gt_kni_temp)
  kr_hb_temp = kr_hb_raw[[i]] %>%
    mutate('time' = rep(i, dim(kr_hb_raw[[i]])[1]))
  kr_hb_exp = rbind(kr_hb_exp, kr_hb_temp)
}
gt_kni_raw_long = gt_kni_exp %>%
  select(c(AP, kni, gt, time)) %>%
  gather(key = gene, value = intensity, kni:gt)

kr_hb_raw_long = kr_hb_exp %>%
  select(c(AP, kr, hb, time)) %>%
  gather(key = gene, value = intensity, kr:hb)

ggplot(gt_kni_raw_long, aes(x = AP, y = intensity, color = gene))+
  geom_point(alpha = 0.8, shape = ".") + 
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression Raw FlyEx Data gt-kni') +
  theme(text = element_text(size = 15))

ggplot(kr_hb_raw_long, aes(x = AP, y = intensity, color = gene))+
  geom_point(alpha = 0.8, shape = ".") + 
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression Raw FlyEx Data kr-hb') +
  theme(text = element_text(size = 15))

# -----------------------------------------------------------------------------------------------------
# 2. clean up data: converting into 100 bins
# -----------------------------------------------------------------------------------------------------
gt_kni_bin = map(gt_kni_raw, Strip_Data, sep = 1, gen1 = "gt", gen2 = "kni")
kr_hb_bin = map(kr_hb_raw, Strip_Data, sep = 1, gen1 = "kr", gen2 = "hb")

# finally assemble the binned data into one complete dataframe
droso_gap_bin = data.frame()
for (i in 1:8){
  gap_temp = cbind('AP' = 4:98, 
                   gt_kni_bin[[i]][4:98, ], kr_hb_bin[[i]][4:98, ], 
                   'time' = rep(i, dim(gt_kni_bin[[i]][4:98, ])[1]))
  droso_gap_bin = rbind(droso_gap_bin, gap_temp)
}

# long df for binned intensity value
droso_gap_bin_long = droso_gap_bin %>%
  gather(key = gene, value = intensity, gt:hb)

# quick plot of the binned intensity results
ggplot(droso_gap_bin_long, aes(x = AP, y = intensity, color = gene)) +
  geom_line(alpha = 0.7) + 
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression Binned FlyEx Data') +
  theme(text = element_text(size = 25))

# -----------------------------------------------------------------------------------------------------
# converting intensity into Boolean Value (1/0 = on/off)
# -----------------------------------------------------------------------------------------------------
gene_mean = map(droso_gap_bin[, 2:5], mean)

droso_Boolean = droso_gap_bin %>%
  mutate(
    'hb_Boolean' = case_when(
      hb > gene_mean['hb'] ~ 1,
      hb < gene_mean['hb'] ~ 0), 
    'kr_Boolean' = case_when(
      kr > gene_mean['kr'] ~ 1,
      kr < gene_mean['kr'] ~ 0),
    'gt_Boolean' = case_when(
      gt > gene_mean['gt'] ~ 1,
      gt < gene_mean['gt'] ~ 0),
    'kni_Boolean' = case_when(
      kni > gene_mean['kni'] ~ 1,
      kni < gene_mean['kni'] ~ 0)
  ) %>%
  select(c(AP, time:kni_Boolean))

# -----------------------------------------------------------------------------------------------------
# 3. output cleaned data files into select directory
# -----------------------------------------------------------------------------------------------------
dir_out = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Github/Drosophila_gap_gene/FlyEx_droso_data"
#  Macintosh HD⁩ ▸ ⁨Users⁩ ▸ ⁨jianhongchen⁩ ▸ ⁨Desktop⁩ ▸ ⁨UCSC_AMS_PhD⁩ ▸ ⁨Github⁩ ▸ ⁨Drosophila_gap_gene⁩
write_csv(droso_Boolean, file.path(dir_out, "droso_Boolean.csv"))
write_csv(droso_gap_bin, file.path(dir_out, "droso_Bin.csv"))
  
  
  
  
  
  
  