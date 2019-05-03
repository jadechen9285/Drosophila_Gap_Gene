# -----------------------------------------
# Markov_TP. R
#   build Markov model and utilize transitional probabilities(TP) to characterize the dynamic of 
#   biologcal spacial-temporal data set (drosophila gap gene network), then make prediction as a 
#   data driven model
# -----------------------------------------
# different section of the codes:
#           Section I:  load libraries and some minor set up.
#           Section II:  Functions.
#           Section III:  Global Variables
#           Section IV:  Actual analysis/method:
#                      1. loading dataset 
#                      2. cleaning up data & analysis
#                      3. simulation 
#                      4. spacial decomposion
#
# -----------------------------------------
# Author:Jianhong Chen
# Date: 04-23-2019
# -----------------------------------------


#########################################################################################
# Section I.  Load the required packages and some minor initial setup :
#########################################################################################

library('ggplot2')
library('tidyverse') # inlcude "tidyr", "readr", "dplyr", "purrr"
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
# Section II.  Created & customized functions specific for this analysis :
#########################################################################################
TP_notation = function(x, dict = state_name) {
  # a function that conver state into 1-16 notation for convienience (built for 'lapply' or 'map' funciton)
  # input:  x: a given state
  #         dict: provided list that acts as Python dictionary object in this case for the given states
  return(which(x == state_name))
}

Spacial_Trans_Prob = function(df){
  # count all the states and calculate the transitional probabilities of the given Boolean state of 
  #   gap genes in selected AP position segments. 
  #       input: df = droso_Boolean (Boolean value of each indiviual gap gene), also need to preprocessed
  #                   based on clustered results (if necessary). 
  #       output: a list of matrices: [[1]] = overall state matrix ; [[2]] = trans_prob matrix ;
  #                                   [[3]] = trans_prob_counter matrix
  
  #1. create the right dataframe for counting the transitional state
  df_state = df %>%
    mutate('gt' = as.character(gt_Boolean), # convert flow data typy into string
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
  row.names(df_state_mat) = unique(df_state$AP) # row = AP position
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
  trans_prob_mat[is.na(trans_prob_mat)] = 0 # replace NA with 0 
  
  return(list(df_state_mat, trans_prob_mat, trans_counter_mat)) # return 3 matrices: state; TP; TP freq.
}

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

Avg_Simulation = function(ge, df, n){
  # average over all time and AP position for specific gene, created for map function mainly
  # input: ge = specific gene(string)
  #        df = dataframe/matrix of simluation results
  #        n = number of simulation
  #        begin = beginning of the AP-position's index
  #        end = end of the AP-position;s index. 
  # output: gene_avg_t = averaged results through all time and AP position(list)
  
  gene_avg_x = gene_avg_t = list()
  for (t in 1:8){
    for (x in unique(df$AP)){
      sim_subset = filter(df, AP == x, time ==t, gene == ge)
      gene_avg_x[[x]] = sum(sim_subset$Boolean)/n #avg at time t through all AP position
    }
    gene_avg_t[[t]] = unlist(gene_avg_x) # avg through all time and all AP
  }
  return(gene_avg_t)
}

Spacial_Sim= function(df, tp, nn){
  
  # perform trans_probability simulation for any given n simulation and AP positions;
  #   then average through n simulation to smooth out all of the randomness.
  # inputs:   df = droso_state_mat; tp = trans_prob
  #           nn = # of simulations; 
  # ouputs:   df_sim_n = a huge matrix/df that contains all the simulation results in genes Boolean
  
  df_n_list  = sim_label = list()
  for (n in 1:nn) {
    df_temp = matrix(0, length(unique(row.names(df))), 8)
    df_temp[, 1] = df[, 1]
    for (t in 1:7){
      df_temp[, t+1] = map_int(df_temp[, t], tp = tp, TP_Calculation)
    }
    # renaming and reorginazing
    row.names(df_temp) = unique(row.names(df))
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
  
  df_avgsim_all = map(gene, Avg_Simulation, df = df_sim_n, n = nn)
  
  # create a time label for data merging during averaging process
  t_label = c()
  for (t in 1:8){
    t_label = c(t_label, rep(t, length(unique(row.names(df)))))
  }
  
  df_avgsim = data.frame('AP' = as.integer(rep(unique(row.names(df)), 8)), 'time' = t_label) %>%
    mutate('gt_avg' = unlist(df_avgsim_all[[1]]),
           'kr_avg' = unlist(df_avgsim_all[[2]]),
           'kni_avg' = unlist(df_avgsim_all[[3]]),
           'hb_avg' = unlist(df_avgsim_all[[4]]) )
  return(df_avgsim)
}
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

# the three functions for customed clustering
Cal_Clusters = function(x, cen){
  # results = map(cen, function(c){norm(as.data.frame(c)-x, type = "2")})
  # cluster = which(min(unlist(results)) == results)
  
  # inputs: x = matrix/df of the clustering data files
  N = dim(x)[1]
  cluster = 1:N
  for (i in 1:N ) {
    results = map_dbl(cen, function(c){norm(c - x[i,], type = "2")})
    cluster[i] = which(min(unlist(results)) == results)
  }
  return(cluster)
}
Update_Centroid = function(k, df,gene){
  df_sub = filter(df, cluster == k)
  cen_new = colMeans(df_sub)
  return(cen_new[gene])
}
Jkmeans = function(df, k){
  # inputs: df = raw input dataframe for the cluster data
  
  # process AP and gene data separately
  
  # generate random centroids
  cen = sample_n(as.data.frame(df), size = k)
  cen_list = list()
  for (i in 1:k){
    cen_list[[i]] = cen[i,]
  }
  
  c_error = rep(10000, k) # cluster error
  p_error = 1e-16 # machine error
  n = 0 # counter
  
  while (norm(c_error, type = "2") > p_error){
    # determine the cluster based on the closest distance
    cluster = Cal_Clusters(df, cen = cen_list)
    df_cl = mutate(as.data.frame(df), "cluster" = cluster)
    
    # update the new centroids by calculating the mean value of each cluster
    cen_new = t(map_dfc(unique(cluster), Update_Centroid, df = df_cl, 
                        gene = names(as.data.frame(df))))
    names(cen_new) = names(as.data.frame(df))
    # calculating the centroids error
    c_error = abs(cen_new - cen)
    # overwrite old centroids with new ones
    cen = cen_new
    n = n + 1 
    print(n)
    print(norm(c_error, type = "2"))
  }
  
  df_cl = mutate(as.data.frame(df), "cluster" = as.factor(cluster))
  
  return(cluster)
}


#########################################################################################
# Section III.  Global Variables
#########################################################################################
gene = c("gt", "kr", "kni", "hb")
cor_index = list(c(2,3),c(2,4),c(2,5),c(3,4),c(3,5),c(4,5)) 
# crossponding label: c(gt-kr, gt-kni, gt-hb, kr-kni, kr-hb, kni-hb)
cor_label = c("gt_kr", 'gt_kni', 'gt_hb', 'kr_kni', 'kr_hb', 'kni_hb')
state_name = c('0 0 0 0', '1 0 0 0', '0 1 0 0', '0 0 1 0', '0 0 0 1',
               '1 1 0 0', '1 0 1 0', '1 0 0 1', '0 1 1 0', '0 1 0 1',
               '0 0 1 1', '1 1 1 0', '1 0 1 1', '1 1 0 1', '0 1 1 1',
               '1 1 1 1' ) # works similiar to Python Dictionary object in this case
nn = 500 # number of simulation

#########################################################################################
# Section IV.  Actual Analysis/Method
#########################################################################################
# -----------------------------------------------------------------------------------------------------
#  1. load the FlyEx data files
# -----------------------------------------------------------------------------------------------------
# current directory that contains the data files:
dir_wd = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Github/Drosophila_gap_gene/FlyEx_droso_data"
setwd(dir_wd)
droso_bin = read_csv("droso_Bin.csv")
droso_Boolean = read_csv("droso_Boolean.csv")
WT_prob = read_csv("droso_WT_prob")
# -----------------------------------------------------------------------------------------------------
#  2. spacial decomposing by correlation based distance clustering
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------
# data preprogressing: obtained corrleation based distance using Spearman Rank method
droso_cluster = droso_bin[, -6]
xx = unique(droso_cluster$AP)

droso_cor = data.frame(row.names = xx)
for (j in 1:length(cor_index)){
  for (i in xx){# i starts with 4 to 98
    df_sub = filter(droso_cluster, AP == i)
    cor_sub=  cor(df_sub[, cor_index[[j]]], method = "spearman")[1,2] # only extra one single corrlation value from the result
    droso_cor[i-3,j] =cor_sub
  }
}
names(droso_cor) = cor_label
droso_cor = mutate(droso_cor, "AP" = xx)
droso_cor_long = gather(droso_cor, key = "gene", value = "Correlation", -AP)

# visulize the distance matrix
ggplot(droso_cor_long, aes(x = AP, y = Correlation)) +
  geom_point(alpha = 0.7) + 
  facet_wrap(.~gene, nrow  = 2) + 
  ylab('gt-kr corrleation') +
  xlab("AP position") +
  ggtitle('Spearman Correlation Distance') +
  theme(text = element_text(size = 25))

# ----------------------------------------------

# optimization k-value of the clustering algorithm 
total_withinSS_km = map_dbl(1:20, function(k){
  model_km = kmeans(droso_cor, centers = k)
  model_km$tot.withinss
})

km_elbow = data.frame( 'k' = 1:20, 'tot_withinSS' = total_withinSS_km)
# the plot implies k = 2 as the optimized value
ggplot(km_elbow, aes(x = k, y = tot_withinSS)) +
  ggtitle("elbow plot for Spearman Correlation Distance") +
  geom_line() +
  scale_x_continuous(breaks = 1:20) + 
  theme(text = element_text(size = 25))

# perform clustering 
cl_exa = kmeans(droso_cor, centers = 5)
cl_jian2 = Jkmeans(droso_cor, k =2)
cl_jian3 = Jkmeans(droso_cor, k =3)
cl_jian4 = Jkmeans(droso_cor, k =4)
cl_jian = Jkmeans(droso_cor, k=5)

df_exa = mutate(filter(droso_bin, time == 3), "cluster" = as.factor(cl_exa$cluster))
df_jian = mutate(filter(droso_bin, time == 3), "cluster" = as.factor(cl_jian))


ggplot(df_jian, aes(x = AP, y = gt, color = cluster )) + 
  geom_point(alpha = 0.6)+
  ylab('Intensity') +
  xlab("AP position") +
  ggtitle('Jian Clustering Method with "Spearman" Correlation Distance') +
  theme(text = element_text(size = 25))

ggplot(df_jian, aes(x = gt, y = kr, color = cluster )) + 
  geom_point(alpha = 0.6) +
  ylab('kr') +
  xlab("gt") +
  ggtitle('Jian Clustering with "Spearman" Correlation Distance') +
  theme(text = element_text(size = 25))

# -----------------------------------------------------------------------------------------------------
#  3. Calculate TP and simulation
# -----------------------------------------------------------------------------------------------------
# calculate the Trans. Probability for clustered data set
jian_Boolean = mutate(droso_Boolean, cluster = rep(cl_jian, 8))
jian_Boolean2 = mutate(droso_Boolean, cluster = rep(cl_jian2, 8))
jian_Boolean3 = mutate(droso_Boolean, cluster = rep(cl_jian3, 8))
jian_Boolean4= mutate(droso_Boolean, cluster = rep(cl_jian4, 8))

jian_tp = list()
for (c in 1:5){
  jian_tp[[c]] = Spacial_Trans_Prob(filter(jian_Boolean, cluster == c))
}

# ----------------------------------------------
# k = 2
jian_tp2 = list()
jian_sim2 = list()

for (c in 1:2){
  jian_tp2[[c]] = Spacial_Trans_Prob(filter(jian_Boolean2, cluster == c))
  jian_sim2[[c]] = Spacial_Sim(df = jian_tp2[[c]][[1]], tp = jian_tp2[[c]][[2]], nn)
}
# assemble all the clustered simulation into one complete results
jian_sim2_final = data.frame()
jian_sim2_final = rbind(filter(jian_sim2[[1]], time == 1), filter(jian_sim2[[2]], time ==1))

for (t in 2:8){
  jian_sim2_final = rbind(jian_sim2_final,
                   filter(jian_sim2[[1]], time == t), filter(jian_sim2[[2]], time ==t))
}
names(jian_sim2_final) = c("AP", "time", "gt_jian2", "kr_jian2", "kni_jian2", "hb_jian2")

# ----------------------------------------------


# ----------------------------------------------
# k = 3
jian_tp3 = list()
jian_sim3 = list()

for (c in 1:3){
  jian_tp3[[c]] = Spacial_Trans_Prob(filter(jian_Boolean3, cluster == c))
  jian_sim3[[c]] = Spacial_Sim(df = jian_tp3[[c]][[1]], tp = jian_tp3[[c]][[2]], nn)
}
# assemble all the clustered simulation into one complete results
jian_sim3_final = data.frame()
jian_sim3_final = rbind(filter(jian_sim3[[3]], time == 1), filter(jian_sim3[[1]], time ==1),
                        filter(jian_sim3[[2]], time ==1))

for (t in 2:8){
  jian_sim3_final = rbind(jian_sim3_final,
                          filter(jian_sim3[[3]], time == t), filter(jian_sim3[[1]], time ==t),
                          filter(jian_sim3[[2]], time == t))
}
names(jian_sim3_final) = c("AP", "time", "gt_jian3", "kr_jian3", "kni_jian3", "hb_jian3")
# ----------------------------------------------


# ----------------------------------------------
# k = 4
jian_tp4 = list()
jian_sim4 = list()

for (c in 1:4){
  jian_tp4[[c]] = Spacial_Trans_Prob(filter(jian_Boolean4, cluster == c))
  jian_sim4[[c]] = Spacial_Sim(df = jian_tp4[[c]][[1]], tp = jian_tp4[[c]][[2]], nn)
}
# assemble all the clustered simulation into one complete results
jian_sim4_final = data.frame()
jian_sim4_final = rbind(filter(jian_sim4[[4]], time == 1), filter(jian_sim4[[2]], time ==1),
                        filter(jian_sim4[[3]], time ==1), filter(jian_sim4[[1]], time == 1))

for (t in 2:8){
  jian_sim4_final = rbind(jian_sim4_final,
                          filter(jian_sim4[[4]], time == t), filter(jian_sim4[[2]], time ==t),
                          filter(jian_sim4[[3]], time == t), filter(jian_sim4[[1]], time ==t))
}
names(jian_sim4_final) = c("AP", "time", "gt_jian4", "kr_jian4", "kni_jian4", "hb_jian4")
# ----------------------------------------------


jian_c1 = Spacial_Trans_Prob(filter(jian_Boolean, cluster == 1))
jian_c2 = Spacial_Trans_Prob(filter(jian_Boolean, cluster == 2))
jian_c3 = Spacial_Trans_Prob(filter(jian_Boolean, cluster == 3))
jian_c4 = Spacial_Trans_Prob(filter(jian_Boolean, cluster == 4))
jian_c5 = Spacial_Trans_Prob(filter(jian_Boolean, cluster == 5))

full_tp = Spacial_Trans_Prob(droso_Boolean)

# perform simulation for each of the cluster 
jian_c1_sim = Spacial_Sim(df = jian_c1[[1]], tp = jian_c1[[2]], nn)
jian_c2_sim = Spacial_Sim(df = jian_c2[[1]], tp = jian_c2[[2]], nn)
jian_c3_sim = Spacial_Sim(df = jian_c3[[1]], tp = jian_c3[[2]], nn)
jian_c4_sim = Spacial_Sim(df = jian_c4[[1]], tp = jian_c4[[2]], nn)
jian_c5_sim = Spacial_Sim(df = jian_c5[[1]], tp = jian_c5[[2]], nn)

full_sim = Spacial_Sim(full_tp[[1]], full_tp[[2]], nn)

# assemble all the clustered simulation into one complete results
jian_sim = data.frame()
jian_sim = rbind(filter(jian_c5_sim, time == 1), filter(jian_c4_sim, time ==1),
                       filter(jian_c2_sim, time == 1), filter(jian_c1_sim, time == 1),
                      filter(jian_c3_sim, time == 1))

for (t in 2:8){
  jian_sim = rbind(jian_sim,
                         filter(jian_c5_sim, time == t), filter(jian_c4_sim, time ==t),
                         filter(jian_c2_sim, time == t), filter(jian_c1_sim, time == t),
                          filter(jian_c3_sim, time == t))
}
names(jian_sim) = c("AP", "time", "gt_jian", "kr_jian", "kni_jian", "hb_jian")

sim_long = gather(jian_sim4_final, key = "gene", value = "avg_prob", gt_jian4:hb_jian4)

ggplot(sim_long, aes(x = AP, y = avg_prob)) +
  geom_line() +
  facet_grid(time ~ gene) +
  ylab('Avg Probabilities') +
  xlab("AP position") +
  ggtitle('Jian Clustering "Spearman" Simulation Results for n = 500') +
  theme(text = element_text(size = 25))

# compare results between different simulations and experimental data
# --------
# construct the dataframe that contain all the simulation results together 
source = c(rep('FullAP_sim', 3040),rep('jian_k2', 3040),rep('jian_k3', 3040),
           rep('jian_k4', 3040), rep('jian_k5', 3040), rep('WT', 3040))
gene_label = rep(c(rep('gt', 760), rep('kr', 760), rep('kni', 760), rep('hb', 760)), 6)
droso_compare_prob = full_sim %>%
  cbind(jian_sim2_final[, names(jian_sim2_final)[3:6]]) %>%
  cbind(jian_sim3_final[, names(jian_sim3_final)[3:6]]) %>%
  cbind(jian_sim4_final[, names(jian_sim4_final)[3:6]]) %>%
  cbind(jian_sim[, names(jian_sim)[3:6]]) %>%
  cbind(WT_prob[, c('gt', 'kr', 'kni', 'hb')]) %>%
  gather(key = 'gene', value = 'value', gt_avg:hb) %>%
  mutate('source' = source, 'gene_label' = gene_label)


# visulized the probabilities comparison results:
ggplot(droso_compare_prob, aes(x = AP, y = value)) +
  geom_line(alpha = 0.8, aes(color = source)) + 
  facet_grid(time~gene_label) + 
  xlab("AP position") +
  ylab("Probability") + 
  ggtitle("Spearman Distance Simulation vs Experimental in Prob N = 500") + 
  theme(text = element_text(size = 22))



