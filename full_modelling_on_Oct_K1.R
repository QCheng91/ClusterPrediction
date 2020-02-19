#Using Oct K1 to do on the cell clustering
#Jan 21 2020
#Including clustering, merging and knn prediction

library(flowCore)
library(tidyverse)
library(stringr)

#Reading Oct K1 data and extract marker counts
get_pro_data <- function() {
  
  #K1
  pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c01_K1_02_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c02_K1_02_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c03_K1_02_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c04_K1_02_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c05_K1_02_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c06_K1_02_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c07_K1_02_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c08_K1_02_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c09_K1_02_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c10_K1_02_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c15_K1_02_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c16_K1_02_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)
  # pro13 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c19_K1_02_0_0_MNCs_Ungated_Ungated.fcs", transformation=FALSE)
  # pro14 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c17_K1_02_0_0_Jurkat_Ungated_Ungated.fcs", transformation=FALSE)
  
  pro1_df = exprs(pro1)
  pro2_df = exprs(pro2)
  pro3_df = exprs(pro3)
  pro4_df = exprs(pro4)
  pro5_df = exprs(pro5)
  pro6_df = exprs(pro6)
  pro7_df = exprs(pro7)
  pro8_df = exprs(pro8)
  pro9_df = exprs(pro9)
  pro10_df = exprs(pro10)
  pro11_df = exprs(pro11)
  pro12_df = exprs(pro12)
  # pro13_df = exprs(pro13)
  # pro14_df = exprs(pro14)
  
  g1 = data.frame(pro1_df)
  g2 = data.frame(pro2_df)
  g3 = data.frame(pro3_df)
  g4 = data.frame(pro4_df)
  g5 = data.frame(pro5_df)
  g6 = data.frame(pro6_df)
  g7 = data.frame(pro7_df)
  g8 = data.frame(pro8_df)
  g9 = data.frame(pro9_df)
  g10 = data.frame(pro10_df)
  g11 = data.frame(pro11_df)
  g12 = data.frame(pro12_df)
  # g13 = data.frame(pro13_df)
  # g14 = data.frame(pro14_df)
  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", 
                 "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4")
  
  other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "MEF2C")
  
  cell_cycle = c("Time", "p_Rb", "cyclin_B1", "IdU")
  
  DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")
  
  barcodes = c("Pd102", "Pd104", "Pd105", "Pd106", "Pd108", "Pd110")
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di", "Sm147Di", "Yb171Di",
              
              "Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di", "Nd146Di",
              
              "Time", "Nd150Di", "Er168Di", "I127Di",
              
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di", 
              
              "Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di") 
  
  pro_marker = c(sur_marker, other_marker, tran_factor, cell_cycle, DNA, barcodes)
  
  seq=c(1:47)
  
  day_levels = paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20))
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    protein1 = bind_rows( g1[marker] %>% mutate(Day = "0", col_id = 1:nrow(g1)), 
                          g2[marker]%>% mutate(Day = "2", col_id = 1:nrow(g2)),
                          g3[marker]%>% mutate(Day = "4", col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(Day = "6", col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(Day = "8", col_id = 1:nrow(g5)),
                          g6[marker]%>% mutate(Day = "10", col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(Day = "11", col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(Day = "12", col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(Day = "14", col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(Day = "16", col_id = 1:nrow(g10)), 
                          g11[marker]%>% mutate(Day = "18", col_id = 1:nrow(g11)),
                          g12[marker]%>% mutate(Day = "20", col_id = 1:nrow(g12))) %>% 
      # g13[marker]%>% mutate(Day = "MNCs", col_id = 1:nrow(g13)),
      # g14[marker]%>% mutate(Day = "Jurkats", col_id = 1:nrow(g14))) %>%
      mutate(Day = factor(Day, levels = day_levels))
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
    
  }  
  
  return(prot_list)
}

all_data = as_tibble(bind_rows(get_pro_data()))

all_data = all_data %>% mutate(cell = factor(paste(cell_id, obstime)))
all_data = all_data %>% mutate(capture = obstime)

gene_data0 = all_data %>% spread(cell.type, Intensity) 

#####gating: filter non-cells######
mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  # reg = lm(data = gated_data, Ir191~Ir193)
  (predict(reg) + ir1)/2
}

gene_data = gene_data0 %>% filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

#check number of cells per day point
cell_count = gene_data %>% group_by(obstime) %>% summarise(count = n())

#extract the interested markers (27+3)
interested_markers = c("CD34", "CD71", "CD45RA", "CD123", "CD49f", 
                       "CD90", "CD44", "CD41", "CD235ab", "CD36", 
                       "CD38", "HBA", "GATA1", "PU1", "ATRX", 
                       "c_Myc", "KLF1", "TAL1", "GATA2", "RUNX1",
                       "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B",
                       "FLI1", "NFE2p45",
                       "mean_Ir", #DNA contents
                       "IdU", "p_Rb") #Cell cycle markers

gated_data = gene_data %>%
  mutate(Day = obstime) %>% 
  select(Day, interested_markers) %>%
  mutate_if(is.numeric, ~ floor(.x)) %>% 
  filter((!Day %in% c("MNCs", "Jurkats")))

######prepare data for clustering model########
prot_name = gated_data %>% select(-Day) %>% colnames() %>% as.character()

prot_list = data.frame(prot_name = prot_name, prot_index = 1:length(prot_name)) %>% as.tibble() %>% 
  mutate(prot_name = as.character(prot_name))

#summarise the gated_data
gated_data %>% group_by(Day) %>% summarise(count=n())

#sub-smaple the cells
cells_per_day = 1000
abv_max = gated_data %>% group_by(Day) %>% filter(n()<=cells_per_day) %>% ungroup()
sub_data = gated_data %>% group_by(Day) %>% filter(n()>cells_per_day) %>% sample_n(cells_per_day) %>% ungroup() %>%
  bind_rows(abv_max) %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  arrange(Day) 

x_data   = sub_data %>% select(-Day) %>% as.matrix()
vec_d    = sub_data %>% mutate(new_indx = as.numeric(factor(as.character(Day), levels = unique(as.numeric(Day))))) %>% pull(new_indx)
indx_tp  = sub_data %>% group_by(Day) %>% summarise(num_cells = n()) %>% mutate(num_cells  = cumsum(num_cells)) %>% 
  mutate(from = c(0, num_cells[-n()]) + 1) %>% select(from, num_cells) %>% as.matrix()

x_indx = x_data %>% as_tibble() %>% mutate_all(~as.numeric(as.factor(.x)))
x_unique_aux = x_data %>% as_tibble() %>% map2(1:ncol(.), ~ as.factor(.x) %>% levels() %>% 
                                                 tibble(unique_count = ., prot =.y)) %>% bind_rows()
unique_prot_x = x_unique_aux %>% mutate(unique_count = as.numeric(unique_count)) %>% pull(unique_count)
indx_prot  = x_unique_aux %>% group_by(prot) %>% summarise(num_vals = n()) %>% mutate(num_vals  = cumsum(num_vals)) %>% 
  mutate(from = c(0, num_vals[-n()]) + 1) %>% select(from, num_vals) %>% as.matrix()

#define the data list
data_cyto = list(c = 30, k = ncol(x_data), N = nrow(x_data), x = x_data, d = max(vec_d), 
                 day_index = vec_d, indx_tp = indx_tp,
                 num_unique_counts = length(unique_prot_x),
                 unique_prot_x = unique_prot_x,
                 indx_prot = indx_prot,
                 x_indx = x_indx)

vec_d %>% unique()

save(sub_data, file = "Results_by_date/Results2020/sampled_sub_data_oct_k1.rda")
#################stan model for clustering, run in parallel############
library(rstan)

options(mc.cores = parallel::detectCores())
doMC::registerDoMC(cores = 14)
rstan_options(auto_write = TRUE)

model   = stan_model("~/Software/cluster_model/mix_nb_model_vec_V10.stan") #mix_nb_model_vec_V6.stan

#vb a function in stan model
clusters = c(30)
runs = c(1)

library(doParallel)

foreach(r = runs)%dopar%{
  
  for(number in clusters) {

    prior = kmeans(asinh(x_data), number, iter.max = 1000)
    
    pri_var = x_data %>% as.tibble() %>% mutate(cluster = prior$cluster) %>% 
      gather(prot, vals, -cluster) %>% group_by(prot, cluster) %>% 
      summarise(variance = var(vals)) %>% spread(prot, variance) %>% 
      select(prot_name) %>% as.matrix()
    
    mu_s_cp = prior$centers %>% sinh() %>% as.tibble() %>% 
      mutate_all(~ replace(., .<0.1, 0.11)) %>% 
      as.matrix(ncol = length(prot_name))
    
    phi_s_cp = ((mu_s_cp * mu_s_cp) / abs(pri_var - mu_s_cp)) %>%
      as.data.frame() %>% #/ abs(pri_var - mu_s_cp) 
      replace(.<=0.1, 0.11) %>% 
      as.matrix(ncol = length(prot_name))
    
    init_vals0 = list(mu_s_cp = mu_s_cp, phi_s_cp = phi_s_cp)
    
    data_cyto = list(c = number, k = ncol(x_data), N = nrow(x_data), x = x_data, d = max(vec_d), 
                     day_index = vec_d, indx_tp = indx_tp,
                     num_unique_counts = length(unique_prot_x),
                     unique_prot_x = unique_prot_x,
                     indx_prot = indx_prot,
                     x_indx = x_indx)
    
    start_time <- Sys.time()
    map_est4   = optimizing(model, data = data_cyto, iter = 2000, as_vector = F, init = init_vals0)
    #map_est4   = vb(model, data = data_cyto, algorithm = "fullrank")
    #print(start_time - Sys.time())
    
    direct = paste("Results_by_date/Results2020/", "full_protein_Ir_Feb2020", number, "_", r, ".rda", sep="")
    save(map_est4, file = direct)
    
  }
  
  return(start_time - Sys.time())
}

############Post-processing the result#########################
load(file = "Results_by_date/Dec2019/sampled_sub_data_oct_k1.rda")
load(file = "Results_by_date/Results2020/full_protein_Ir_Feb202030_1.rda") #map_est4 #full_protein_Ir_init_c_test25_1.rda Dec2019/full_protein_Ir_init_c_test25_1.rda
load(file = "Results_by_date/Dec2019/full_protein_Ir_init_c_test25_1.rda")
run_fit = map_est4

day_vals = sub_data$Day %>% unique()
day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)

c = 25

vals = run_fit$par 

time = tibble(Day = day_index_0$day_index,  real_time= day_index_0$Day) 

cluster = vals$theta %>% t() %>% as.tibble() %>%
  rename(vals = V1) %>%
  mutate(cluster = 1:c)

phi = vals$phi_s_cp %>% as.tibble() %>% 
  mutate(cluster = 1:c) %>% 
  gather(protein, std, -cluster) %>%
  mutate(protein = as.double(str_remove_all(protein, "V")))

signal = vals$mu_s_cp %>% as.tibble() %>% 
  mutate(cluster = 1:c) %>% 
  gather(protein, vals, -cluster) %>%
  mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
  left_join(prot_list, by = c("protein" = "prot_index")) %>% 
  left_join(phi) 

signal = signal %>% 
  group_by(prot_name) %>% 
  mutate(max_vals = max(vals), min_vals = min(vals), mid_vals = mean(vals), sd_vals = sd(vals)) %>% 
  ungroup()

day_vals = gated_data$Day %>% unique()
day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)

mu_c_fit = vals$mu_s_cp %>% t()

#theta0_c_fit = vals$theta0_s_cp %>% t()

phi_c_fit = vals$phi_s_cp %>% t()

theta_fit = vals$theta %>% t()

data_x = gene_data0 %>% 
  filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
  mutate(Day = obstime) %>%
  as.tibble() %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 

data_x1 = data_x[-650918,] %>% select(prot_name) %>% mutate_if(is.numeric, ~ floor(.x)) %>% as.matrix()

day_indx = data_x[-650918,] %>%  pull(day_index)

#calculate_log_likeli.cpp
#####predict cell cluster classifications
all_cells = predict_cluster(data_x1, theta_fit, mu_c_fit, phi_c_fit, day_indx)

#####calculate the probabilities in each cluster
all_cells_prob = predict_cluster_prob(data_x1, theta_fit, mu_c_fit, phi_c_fit, day_indx)

# check any wrong readings 
# data_x[650918,] %>% View()
# all_cells_prob[650918,]

####################plotting the clustering result########################
cluster_prob = all_cells %>% as.tibble()

full_clusters = all_cells %>% as.tibble() %>% unique() %>% pull(V1)

p1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>% 
  #filter(cluster %in% full_clusters) %>% 
  #filter(prot_name %in% key_proteins) %>% 
  ggplot(aes(x = prot_name, y= vals)) + 
  geom_point(aes(colour = prot_name)) +
  scale_y_log10()+
  geom_text(aes(label = prot_name), size=3.0)+
  facet_wrap(~cluster, nrow = 2) +
  theme(legend.position = "none")

p2 = data_x[-650918,] %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(cluster, prob, fill = 0) %>% 
  gather(cluster, prob, -Day) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  ggplot(aes(x = Day, y = prob))+
  geom_point()+
  facet_wrap(~cluster, nrow=2)

gridExtra::grid.arrange(p1, p2, ncol =1)

