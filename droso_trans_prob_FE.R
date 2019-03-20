# -----------------------------------------
# droso_trans_prob_FE. R
#   implement transitional probabilities(Markov-like) methods to analysis 
#   drosophilia gag gene network
#   database: FlyEx
# -----------------------------------------
# different section of the codes:
#           I. setup and loading dataset 
#           II. cleaning up data & analysis
#           III. simulation 
#           IV. spacial decomposion
# -----------------------------------------
# Author:Jianhong Chen
# Date: Jan 20th 2019
# -----------------------------------------


# load required packages

library('ggplot2')
library('tidyr')
library('dplyr')
library('tidyverse')
library('purrr')
library('reshape2') 
library('plotly') # 3D plotting
library("kernlab") # kernel-based machine learning techniques 
library("factoextra") # multivariable analysis in R [hkmean]
library("dendextend") # create dendrogram object for hclustering plot 
library("cluster") # clustering functions
library("fpc") # Density-Based Spatial Clustering and application with Noise
library("dbscan")# Density-Based Spatial Clustering and application with Noise

# -----------------------------------------------------------------------------------------------------
# intial set up
ls() # list all objects 
gc() # force R to release memory of deleted variables 
rm(list = ls()) # delete all variables in the environment
# -----------------------------------------------------------------------------------------------------


#########################################################################################
# I.  setup and loading dataset :
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# load the FlyEx data files
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
  #gt_kni_exp = rbind(gt_kni_exp, gt_kni_temp)
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


#########################################################################################
# II.  cleanup & analysis :
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# clean up data: converting into 100 bins
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
  theme(text = element_text(size = 15))

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

droso_Boolean_long = droso_Boolean %>%
  gather(key = gene, value = Boolean, hb_Boolean:kni_Boolean)
  
# quick plot of Boolean Value 
ggplot(droso_Boolean_long, aes(x= AP, y = Boolean, color = gene)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression in Boolean') +
  theme(text = element_text(size = 15))

# -----------------------------------------------------------------------------------------------------
# create gene state: (gt, kr, kni, hb) 
#   calculate the transitional probabilitie between each different state
# Notice: Markov-like behaviro: only the previous timeframe can affects the next one
# -----------------------------------------------------------------------------------------------------
gene = c("gt", "kr", "kni", "hb")
state_name = c('0 0 0 0', '1 0 0 0', '0 1 0 0', '0 0 1 0', '0 0 0 1',
               '1 1 0 0', '1 0 1 0', '1 0 0 1', '0 1 1 0', '0 1 0 1',
               '0 0 1 1', '1 1 1 0', '1 0 1 1', '1 1 0 1', '0 1 1 1',
               '1 1 1 1' ) # works similiar to Python Dictionary object in this case

TP_notation = function(x, dict = state_name) {
  # a function that conver state into 1-16 notation for convienience (built for 'lapply' or 'map' funciton)
  # input:  x: a given state
  #         dict: provided list that acts as Python dictionary object in this case for the given states
  return(which(x == state_name))
}


Spacial_Trans_Prob = function(df, begin, end){
  # counter all the states and calculate the transitional probabilities of the given Boolean state of 
  #   gap genes in selected AP position segments. 
  #       input: df = droso_Boolean (Boolean value of each indiviual gap gene); 
  #              begin = AP starting point; end = AP ending point;
  #       output: a list of matrices: [[1]] = overall state matrix ; [[2]] = trans_prob matrix ;
  #                                   [[3]] = trans_prob_counter matrix
  
  #1. create the right dataframe for counting the transitional state
  df_state = filter(df, (AP >= begin & AP <= end))%>%
    mutate('gt' = as.character(gt_Boolean),
           'kr' = as.character(kr_Boolean),
           'kni' = as.character(kni_Boolean),
           'hb' = as.character(hb_Boolean)) %>%
    mutate('state' = paste(gt, kr, kni, hb)) %>%
    select(c("AP", "time", "state")) 
  
  df_state_int = data_frame("state" = map_int(matrix(unlist(df_state[["state"]])), TP_notation)) %>%
    mutate("time" = df_state[["time"]])
  
  # make droso_state_string into a matrix format
  df_state_mat = subset(df_state_int, time == 1)
  for (t in 2:8){
    df_state_mat = cbind(df_state_mat, subset(df_state_int, time == t))
  }
  df_state_mat = df_state_mat[, c(1,3,5,7,9,11,13,15)] # filter out the state variables only
  row.names(df_state_mat) = begin:end # row = AP position
  colnames(df_state_mat) = 1:8 # col = timeframe
  
  # 2. count the number of the transition
  trans_counter_mat = matrix(0, 16, 16) # row = initial_state, col = trans_state
  for (x in 1:dim(df_state_mat)[1]){
    for (t in 1:7){
      initial_state = df_state_mat[x,t]
      trans_state = df_state_mat[x, t+1]
      trans_counter_mat[initial_state, trans_state] = trans_counter_mat[initial_state, trans_state] + 1
    }
  }
  
  # 3. convert the results into probabilities
  trans_prob_list = list()
  for (i in 1:16){
    trans_prob_temp = map(trans_counter_mat[i, ], function(x){x/sum(trans_counter_mat[i, ])})
    trans_prob_list= append(trans_prob_list, trans_prob_temp)
  }
  trans_prob_mat = matrix(unlist(trans_prob_list), nrow = 16) # row = transitional state, col = initial state
  trans_prob_mat[is.na(trans_prob_mat)] = 0 # replace NA with 0 entries
  
  return(list(df_state_mat, trans_prob_mat, trans_counter_mat))
}

# kernel k-mean clustered result
cluster1_1_trans_state = Spacial_Trans_Prob(droso_Boolean, 4, 23)[[1]]
cluster2_1_trans_state = Spacial_Trans_Prob(droso_Boolean, 24, 39)[[1]]
cluster1_2_trans_state = Spacial_Trans_Prob(droso_Boolean, 40, 46)[[1]]
cluster2_2_trans_state = Spacial_Trans_Prob(droso_Boolean, 47, 69)[[1]]
cluster3_trans_state = Spacial_Trans_Prob(droso_Boolean, 70, 98)[[1]]

cluster1_1_tp = Spacial_Trans_Prob(droso_Boolean, 4, 23)[[2]]
cluster2_1_tp = Spacial_Trans_Prob(droso_Boolean, 24, 39)[[2]]
cluster1_2_tp = Spacial_Trans_Prob(droso_Boolean, 40, 46)[[2]]
cluster2_2_tp = Spacial_Trans_Prob(droso_Boolean, 47, 69)[[2]]
cluster3_tp = Spacial_Trans_Prob(droso_Boolean, 70, 98)[[2]]

cluster1_1_tp_counter = Spacial_Trans_Prob(droso_Boolean, 4, 23)[[3]]
cluster2_1_tp_counter = Spacial_Trans_Prob(droso_Boolean, 24, 39)[[3]]
cluster1_2_tp_counter = Spacial_Trans_Prob(droso_Boolean, 40, 46)[[3]]
cluster2_2_tp_counter = Spacial_Trans_Prob(droso_Boolean, 47, 69)[[3]]
cluster3_tp_counter = Spacial_Trans_Prob(droso_Boolean, 70, 98)[[3]]



full_trans_state = Spacial_Trans_Prob(droso_Boolean, 4,98)[[1]]
full_tp = Spacial_Trans_Prob(droso_Boolean, 4, 98)[[2]]
full_tp_counter = Spacial_Trans_Prob(droso_Boolean, 4, 98)[[3]]

# visulizating the reulst
tp_long = melt(cluster3_tp)

ggplot(tp_long, aes(x = Var1, y = value)) +
  geom_col() +
  facet_wrap(~Var2) +
  ylab('total count') +
  xlab("Transitonal States") +
  ggtitle('Full Range Transitional Probabilities Histogram') +
  theme(text = element_text(size = 20))

ggplot(tp_long, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low="grey90", high="red") +
  xlab('Initial States') +
  ylab("Transitonal States") +
  ggtitle('Cluster 3 Transitional Probabilities Matrix') +
  theme(text = element_text(size = 20))

ggplot(tp_long, aes(x = Var1, y = value)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Var2) +
  ylab('Transitional Probabilities') +
  xlab("Transitonal States") +
  ggtitle('Cluster 3 Transitional Probabilities Barplot') +
  theme(text = element_text(size = 20))

#########################################################################################
# III.  Simulation:
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# Simulate the gap gene interactions with the calculated transitional probabilities
# -----------------------------------------------------------------------------------------------------
TP_Calculation = function(s, tp){
  # Scan through all states and determine the next state based on the calculated transitional probabilities
  #       input: tp = trans_prob; s = initial state
  #       output: next_state = state after applying transitional probabilities
  
  # obtain a cummulative version of trans.prob. matrix for simulation purpose
  trans_prob_cummulative = tp # cummlative version of the trans_prob for simulation purpose
  for (init in 1:16){
    ind_zero = which(tp[, init] == 0) # record the index where the entry ==0
    for (tran in 1:15){
      trans_prob_cummulative[tran+1, init] = trans_prob_cummulative[tran, init] + trans_prob_cummulative[tran+1, init]
    }
    trans_prob_cummulative[ind_zero, init] = 0 # reset the original entries back to zero
  }
  # now scan through each compact state of the given file 
  tp_one = trans_prob_cummulative[,s]
  for (i in tp_one){
    rd = runif(1) # generate random number
    if (rd < i){
      next_state = which(tp_one == i)
      break # it is important to break the for loop here to terminate scanning which change i 
    }
  }
  return(next_state)
} 

Avg_Simulation = function(ge, df, n, begin, end){
  # average over all time and AP position for specific gene, created for map function mainly
  # input: ge = specific gene(string)
  #        df = dataframe/matrix of simluation results
  #        n = number of simulation
  #        begin = beginning of the AP-position's index
  #        ened = end of the AP-position;s index. 
  # output: gene_avg_t = averaged results through all time and AP position(list)
  
  gene_avg_x = gene_avg_t = list()
  for (t in 1:8){
    for (x in begin:end){
      sim_subset = filter(df, AP == x, time ==t, gene == ge)
      gene_avg_x[[x]] = sum(sim_subset$Boolean)/n #avg at time t through all AP position
    }
    gene_avg_t[[t]] = unlist(gene_avg_x) # avg through all time and all AP
  }
  return(gene_avg_t)
}

nn = 100 # determine number of simulation

Spacial_Sim= function(df, tp, nn, begin, end){
  
  # perform trans_probability simulation for any given n simulation and AP positions;
  #   then average through n simulation to smooth out all of the randomness.
  # inputs:   df = droso_state_mat; tp = trans_prob
  #           nn = # of simulations; begin starting AP position; end = ending AP position
  # ouputs:   df_sim_n = a huge matrix/df that contains all the simulation results in genes Boolean
  
  df_n_list  = sim_label = list()
  for (n in 1:nn) {
    df_temp = matrix(0, (end - begin+1), 8) # length(end-begin+1) = length(begin:end)
    df_temp[, 1] = df[, 1]
    for (t in 1:7){
      df_temp[, t+1] = map_int(df_temp[, t], tp = tp, TP_Calculation)
    }
    # renaming and reorginazing
    row.names(df_temp) = begin:end
    df_temp_long = melt(df_temp)
    # decoding the state back to (gt, kr, kni, hb):
    df_de = data_frame('AP' = df_temp_long$Var1, 
                       'time' = df_temp_long$Var2, 
                       'state' = state_name[df_temp_long$value])
    state_split = strsplit(df_de$state, split =' ')
    gt = kr = kni = hb = list()
    for (i in state_split){
      gt = unlist(append(gt, as.integer(i[1])))
      kr = unlist(append(kr, as.integer(i[2])))
      kni = unlist(append(kni, as.integer(i[3])))
      hb = unlist(append(hb, as.integer(i[4])))
    }
    df_final = data.frame(df_de) %>%
      mutate(gt, kr, kni, hb) %>% 
      gather(key = 'gene', value = 'Boolean', gt:hb)
    # store each simulation into a large list
    df_n_list[[n]] = df_final
    # creating 'sim' vector for naming purpose
    sim_label[[n]] = rep(n, dim(df_final)[1])
  }
  
  # now concatenate all simulations into one huge matrix/df:
  df_sim_n = df_n_list[[1]]
  for (n in 2:nn){ #start with 2, because the 1st one is used for initiaion of the matrix, 'droso_sim_n'
    df_sim_n = rbind(df_sim_n, df_n_list[[n]])
  }
  df_sim_n = mutate(df_sim_n, 'sim' = unlist(sim_label))
  
  df_avgsim_all = map(gene, Avg_Simulation, df = df_sim_n, n = nn, begin, end)
  
  # create a time label for data merging during averaging process
  t_label = c()
  for (t in 1:8){
    t_label = c(t_label, rep(t, (end - begin + 1)))
  }
  
  df_avgsim = data.frame('AP' = rep(begin:end, 8), 'time' = t_label) %>%
    mutate('gt_avg' = unlist(df_avgsim_all[[1]]),
           'kr_avg' = unlist(df_avgsim_all[[2]]),
           'kni_avg' = unlist(df_avgsim_all[[3]]),
           'hb_avg' = unlist(df_avgsim_all[[4]]) )
  return(df_avgsim)
}

droso_sim_100 = Spacial_Sim(full_trans_state,tp = full_tp, nn, 4, 98)

c1_1_sim = Spacial_Sim(cluster1_1_trans_state, tp = cluster1_1_tp, nn, 4, 23)
c2_1_sim = Spacial_Sim(cluster2_1_trans_state, tp = cluster2_1_tp, nn, 24, 39)
c1_2_sim = Spacial_Sim(cluster1_2_trans_state, tp = cluster1_2_tp, nn, 40, 46)
c2_2_sim = Spacial_Sim(cluster2_2_trans_state, tp = cluster2_2_tp, nn, 47, 69)
c3_sim = Spacial_Sim(cluster3_trans_state, tp = cluster3_tp, nn, 70, 98)

# stich all the clusters results into one complete one
all_cluster = data.frame(rbind(filter(c1_1_sim, time == 1), filter(c2_1_sim, time == 1),
                         filter(c1_2_sim, time == 1),filter(c2_2_sim, time == 1), filter(c3_sim, time == 1)))
for (t in 2:8){
  all_cluster = rbind(all_cluster,(rbind(filter(c1_1_sim, time == t), filter(c2_1_sim, time == t),
                      filter(c1_2_sim, time == t),filter(c2_2_sim, time == t), filter(c3_sim, time == t))))
}
colnames(all_cluster) = c("AP", 'time', 'gt_avg_cl', 'kr_avg_cl', 'kni_avg_cl', 'hb_avg_cl') # relabeling

# visulization the simulation results:
sim_long = gather(all_cluster, key = 'gene', value = 'avg_prob', gt_avg_cl:hb_avg_cl)

ggplot(sim_long, aes(x = AP, y = avg_prob)) +
  geom_line() +
  facet_grid(time ~ gene) +
  ylab('Avg Probabilities') +
  xlab("AP position") +
  ggtitle('3rd Clustering Simulation Results for n = 100') +
  theme(text = element_text(size = 20))

#-------------------------------------------------------##########

# calculate the WT probabilities: intensity >> Boolean >> probabilities >> bin to 100
Probability_Conversion = function(x, sep,  gene1, gene2){
  # Convert the raw intensity measurements of the gap genes into 
  #     probabilities by counting the Boolean value then average over a strip:
  #   inputs: x = data file;
  #           sep = width of the strip
  #           gen1, gen2 = interested genes
  #   output: a dataframe of converted probabilities of the gap gene
  one_prob = data.frame()
  for (n in seq(0, (100-sep),sep)){
    one = subset(x, (x$AP > n) & (x$AP < (n+sep)))
    gene1_bool = map(one[[gene1]], function(x){case_when(x > gene_mean[gene1]~1, x < gene_mean[gene1] ~ 0)})
    gene2_bool = map(one[[gene2]], function(x){case_when(x > gene_mean[gene2]~1, x < gene_mean[gene2] ~ 0)})
    gene1_prob = sum(unlist(gene1_bool))/length(gene1_bool)
    gene2_prob = sum(unlist(gene2_bool))/length(gene2_bool)
    one_prob = rbind(one_prob, c(gene1_prob, gene2_prob))
  }
  colnames(one_prob) = c(gene1, gene2)
  return(one_prob)
}
gt_kni_prob = map(gt_kni_raw, Probability_Conversion, sep = 1, gene1 = 'gt', gene2 = 'kni')
kr_hb_prob = map(kr_hb_raw, Probability_Conversion, sep = 1, gene1 = 'kr', gene2 = 'hb')

# finally assemble the prob data into one complete dataframe
WT_prob = data.frame()
for (i in 1:8){
  prob_temp = cbind('AP' = 4:98, # only consider data in this range because of 'NAN'
                   gt_kni_prob[[i]][4:98, ], kr_hb_prob[[i]][4:98, ], 
                   'time' = rep(i, dim(gt_kni_prob[[i]][4:98, ])[1]))
  WT_prob = rbind(WT_prob, prob_temp)
}

# Combines normal(whole space sim), clustering sim, and WT results into one complete df for visilization
source = c(rep('sim', 3040),rep('cl_sim', 3040),rep('WT', 3040))
gene_label = rep(c(rep('gt', 760), rep('kr', 760), rep('kni', 760), rep('hb', 760)), 3)
droso_compare_prob = droso_sim_100 %>%
  cbind(all_cluster[, c('gt_avg_cl', 'kr_avg_cl', 'kni_avg_cl', 'hb_avg_cl')]) %>%
  cbind(WT_prob[, c('gt', 'kr', 'kni', 'hb')]) %>%
  gather(key = 'gene', value = 'value', gt_avg:hb) %>%
  mutate('source' = source, 'gene_label' = gene_label)

# visulized the probabilities comparison results:
ggplot(droso_compare_prob, aes(x = AP, y = value)) +
  geom_line(alpha = 0.8, aes(linetype = source)) + 
  facet_grid(time~gene_label) + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid")) +
  xlab("AP position") +
  ylab("Probability") + 
  ggtitle("Simulation vs Experimental in Prob N = 100") + 
  theme(text = element_text(size = 20))

# -----------------------------------------------------------------------------------------------------
# simulation with diffusion D(x_i-1 + x_i + x_j+1) : diffusion >>>>> Trans. Prob. (for every timestep)
# -----------------------------------------------------------------------------------------------------

 
#########################################################################################
# IV.  Spacial Decomposition :
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# Experimental EDA spcaially 
# -----------------------------------------------------------------------------------------------------

# first compare the effect of binning, see if binning change the dynamic of the interaction
ggplot(droso_gap_bin, aes(x = gt, y = kni, color = AP)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~time, nrow =2) +
  ggtitle("gt vs kni in time binned data set ") + 
  theme(text = element_text(size = 15))

ggplot(gt_kni_exp, aes(x = gt, y = kni, color = AP)) +
  geom_point(alpha = 0.8, shape =  '.') +
  facet_wrap(~time, nrow = 2) +
  ggtitle("gt vs kni in time raw data set ") + 
  theme(text = element_text(size = 15))

# binned strong mutual repression pairs interaction 
ggplot(droso_gap_bin, aes(x = gt, y = kr, color = AP)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~time, nrow =2) +
  ggtitle("gt vs kr in time binned data set ") + 
  theme(text = element_text(size = 15))

# create an extra label for 'droso_gap_bin' df to create manual segmentation 
AP_label = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100")
seg_label= c()
for (x in seq(0, 90, 10)){
  seg_temp = filter(droso_gap_bin, time == 1 & AP >= x & AP < x +10)
  seg_label = c(seg_label, rep(AP_label[x/10 + 1], dim(seg_temp)[1]))
}
droso_gap_bin = mutate(droso_gap_bin, "manual_seg" = rep(seg_label, 8))

ggplot(droso_gap_bin, aes(x = gt, y = kr)) + 
  geom_point() + 
  facet_grid(manual_seg ~ time) +
  ggtitle("gt vs kr in time with manual segmentation binned data set ") + 
  theme(text = element_text(size = 15))

ggplot(droso_gap_bin, aes(x = kni, y = hb)) + 
  geom_point() + 
  facet_grid(manual_seg ~ time) +
  ggtitle("kni vs hb in time with manual segmentation binned data set ") + 
  theme(text = element_text(size = 15))

# -----------------------------------------------------------------------------------------------------
# Implement Clustering Analysis Methods: k-mean , kernel k-mean(nonlinear), Hierarchical k-mean (hkmean)
# -----------------------------------------------------------------------------------------------------
# output "droso_gap_bin" into a csv file
# write_csv(droso_gap_bin, path = "gap_gene_intensity")

km_df = scale(droso_gap_bin[ , c(-6, -7)]) # scale the df for clustering purpose

set.seed(78)
# first determine the optimize k-value for the data based on elbow and silhouette method
# use elbow on regular k-mean
total_withinSS_km = map_dbl(1:20, function(k){
  model_km = kmeans(km_df, centers = k)
  model_km$tot.withinss
})

km_elbow = data.frame( 'k' = 1:20, 'tot_withinSS' = total_withinSS_km)
# the plot implies k = 2 as the optimized value
ggplot(km_elbow, aes(x = k, y = tot_withinSS)) +
  ggtitle("elbow plot for gene intensity") +
  geom_line() +
  scale_x_continuous(breaks = 1:20) + 
  theme(text = element_text(size = 20))

# silhouette method for kmean
sil_width = map_dbl(2:20, function(k){
  pam_km = pam(km_df, k = k)
  pam_km$silinfo$avg.width
})
sil_df = data.frame('k' = 2:20, 'avg_sil_width' = sil_width)

# silhouette results : k = 5
ggplot(sil_df, aes(x = k, y = avg_sil_width)) +
  ggtitle("silhouette plot for gene intensity") +
  geom_line() +
  scale_x_continuous(breaks = 2:20) +
  theme(text = element_text(size = 20))
  

# implement k-mean
model_km = kmeans(km_df, centers =5)
droso_km = mutate(droso_gap_bin, "cluster" = model_km$cluster)
                          
# --------------------------------------
# implement kernel k-mean
model_kkm = kkmeans(as.matrix(km_df), centers = 3) # kkmeans only takes in matrix/list data type
droso_kkm = mutate(droso_gap_bin, "cluster" = model_kkm)

# implement Hierachical k-mean: an optimized version for k-mean 
model_hkm = hkmeans(km_df, k = 3)
droso_hkm = mutate(droso_gap_bin, "cluster" = model_hkm$cluster)

# implement Fuzzy clustering: fuzzy c-mean (FCM)
model_fcm = fanny(km_df, k = 3, metric = "euclidean", stand = FALSE)
droso_fcm = mutate(droso_gap_bin, "cluster" = model_fcm$clustering)

# implement density-based spatial clustering and application with noise(DBSCAN)
# 1. need to first find the optimized epsilon-value and minimum points
dbscan::kNNdistplot(km_df, k = 5)
abline(h = 0.5, lty = 2)
# >> k-NN results implies that mini.pts = 5 and epsilon-value = 0.5

# 2. build the density-based cluster model 
model_db = fpc::dbscan(km_df, eps = 0.15, MinPts = 5)


# visualize the dendrogram/tree 
fviz_dend(model_hkm)
fviz_cluster(model_fcm)
fviz_cluster(model_db, data = km_df)

#3D plot for the clustered data through all time.
km_3D <- plot_ly(droso_kkm3, x = ~AP, y = ~gt, z = ~kr,
             marker = list(size = 2, color = ~as.factor(cluster), showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'AP'),
                      yaxis = list(title = 'gt intensity'),
                      zaxis = list(title = 'kr intensity')),
         annotations = list(
           x = 1.15,
           y = 1.05,
           text = 'cluster',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))

# 2D plot of the cluster results
ggplot(droso_fcm, aes(x = AP , y = gt, color = factor(cluster))) + 
  geom_point(alpha =0.6) + 
  xlab("AP") +
  ylab("gt") + 
  ggtitle("Fuzzy c-means AP vs gt Intensity") +
  theme(text = element_text(size = 20))

ggplot(droso_fcm, aes(x = gt , y = kr, color = factor(cluster))) + 
  geom_point(alpha =0.6) + 
  facet_wrap(~time, nrow =2 ) + 
  xlab("gt") +
  ylab("kr") + 
  ggtitle("Fuzzy c-means gt vs kr Intensity") +
  theme(text = element_text(size = 20))

# -------------------------------
# perform simulation with kernel kmean result at time 4 
# -------------------------------
kkm3_t4_c1 = filter(droso_kkm3, time == 4 & cluster == 1)
kkm3_t4_c2 = filter(droso_kkm3, time == 4 & cluster == 2)
kkm3_t4_c3 = filter(droso_kkm3, time == 4 & cluster == 3)


#########################################################################################
# V.  Time Delay :
#########################################################################################
time_label = c('1-2', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8')
gene_lag  = data.frame()
gene = c('gt', 'kr', 'kni', 'hb')
for (g in gene){
  for (t in 1:7){
    gene_t1 = filter(droso_gap_bin, time == t)[g]
    gene_t2 = filter(droso_gap_bin, time == t+1)[g]
    gene_temp = data.frame(gene_t1, gene_t2, rep(time_label[t], 95), rep(g, 95))
    colnames(gene_temp) = c('time1', 'time2', 'time_label', 'gene') 
    gene_lag = rbind(gene_lag, gene_temp)
  }
}

ggplot(gene_lag, aes(x = time1, y = time2))+
  geom_point(alpha = 0.8) + 
  facet_grid(time_label~gene) +
  xlab("Time (t)") +
  ylab("Time (t + 1)") +
  ggtitle("Time Delayed(lag by 1)") +
  theme(text = element_text(size = 15))

 
