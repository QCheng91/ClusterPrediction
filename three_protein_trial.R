#7/8/2019
#new trial for the proof of principal
#3 protein model

library(flowCore)
library(tidyverse)
library(stringr)

#Oct data
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
  # 
  # # #K2
  # pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c01_K2_01_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c02_K2_01_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c03_K2_01_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c04_K2_01_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c05_K2_01_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c06_K2_01_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c07_K2_01_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  # pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c08_K2_01_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c09_K2_01_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c10_K2_01_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  # pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c15_K2_01_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  # pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c20_K2_01_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)
  # 
  #tibble(desc = as.character(pro13@parameters@data$desc), name = names(pro13@parameters@data$desc)) %>% separate(desc, c("Isotope", "Protein"), "_") %>% View()
  
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
  
  cell_cycle = c("Time", "p_Rb", "cyclin_B1", "IdU") #p_HH3
  
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
  
  #day_levels = paste0("", c("MNCs", "Jurkats")) #paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20, "MNCs", "Jurkats"))
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

#####gating: filter non-cells

mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  # reg = lm(data = gated_data, Ir191~Ir193)
  (predict(reg) + ir1)/2
}

gene_data = gene_data0 %>% filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
  filter(mean_Ir >= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5),
         mean_Ir <= sinh(8.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

cell_count = gene_data %>% group_by(obstime) %>% summarise(count = n()) 

cell_count %>%
  ggplot() +
  geom_point(aes(x = obstime, y = count), size = 4, colour = "blue") +
  # geom_errorbar(aes(x = time_point, ymax= exp_total+1.9*exp_std, ymin= exp_total-1.9*exp_std, colour = Exp, group ="Experiment" ), width =0.2)+
  # geom_line(aes(y = `50%`, colour = Exp, linetype = "Modelling"), size = 0.65) + 
  # geom_ribbon(aes(ymin= `2.5%`, ymax =`97.5%`, fill = Exp, size = "Credible interval"), alpha = 0.1) +
  theme_light() +
  theme(text = element_text(size=15), legend.title=element_blank()) +
  #facet_wrap(~Exp)+ #, scales = "free"
  ylab("Total Number of Cells") +
  xlab("Day")

sur_marker = c("CD34", "CD71", "CD45RA", "CD123", "CD49f", 
               "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4",
               "CD36", "CD38")

other_marker = c("HBB", "HBA", "H3")

tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "TAL1", "GATA2",
                "RUNX1", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "Flag_Tag",
                "FLI1", "NFE2p45")

cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")

DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")

interested_markers = c("CD34", "CD71", "CD45RA", "CD123", "CD49f", 
               "CD90", "CD44", "CD41", "CD235ab", "CD36", "CD38", 
               "HBA", "GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "TAL1", 
               "GATA2", "RUNX1", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B",
               "FLI1", "NFE2p45",
               # "Pd102", "Pd104", "Pd105", "Pd106", "Pd108", "Pd110",
               "mean_Ir", #"Ce140", "Pt195", "Pt194",
               "IdU", "p_Rb")

gated_data = gene_data %>%  #express_data0 %>%
  mutate(Day = obstime) %>% 
  select(Day, interested_markers) %>% #-cell_cycle, -HBB, -H3, -DNA, -mean_Ir, -mean_Pt ; -cell_id, -cell, -capture, -obstime
  mutate_if(is.numeric, ~ floor(.x)) %>% 
  filter((!Day %in% c("MNCs", "Jurkats"))) #, IdU, p_Rb

# gated_data = gene_data %>%  #express_data0 %>%
#   mutate(Day = obstime) %>% 
#   select(-cell_id, -cell, -capture, -obstime, -cell_cycle, -HBB, -H3, -DNA, -mean_Ir, -mean_Pt) %>% #
#   mutate_if(is.numeric, ~ floor(.x)) %>% 
#   filter((!Day %in% c("MNCs", "Jurkats")))

prot_name = gated_data %>% select(-Day) %>% colnames() %>% as.character()

prot_list = data.frame(prot_name = prot_name, prot_index = 1:length(prot_name)) %>% as.tibble() %>% 
  mutate(prot_name = as.character(prot_name))

#summarise the gated_data
gated_data %>% group_by(Day) %>% summarise(count=n())

cells_per_day = 2500
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

data_cyto = list(c = 25, k = ncol(x_data), N = nrow(x_data), x = x_data, d = max(vec_d), 
                 day_index = vec_d, indx_tp = indx_tp,
                 num_unique_counts = length(unique_prot_x),
                 unique_prot_x = unique_prot_x,
                 indx_prot = indx_prot,
                 x_indx = x_indx)

vec_d %>% unique()
############################################
library("cytofkit")
library(tidyverse)

tsne_model_1 = Rtsne(as.matrix(express_data), check_duplicates=FALSE, pca=F, perplexity=30, theta=0.5, dims=2)
d_tsne_1 = as.data.frame(tsne_model_1$Y)

ggplot(d_tsne_1, aes(x=V1, y=V2)) +
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")

express_data = sub_data %>% select(-Day) %>% asinh() #%>% as.matrix() gated_data %>% group_by(Day) %>% sample_n(5000) %>% ungroup()
cluster_PhenoGraph = cytof_cluster(xdata = express_data, method = "Rphenograph")
cluster_tsne = cytof_dimReduction(data = express_data, method = "tsne")

gated_data %>% group_by(Day) %>% sample_n(5000) %>% ungroup() %>% select(-Day) %>% asinh()  %>% 
  mutate(tsne_1 = cluster_tsne[,1], tsne_2 = cluster_tsne[,2], pheno_cluster = as.factor(cluster_PhenoGraph)) %>%
  #mutate(CEBPa_flag = if_else(CEBPa >= 4.0, 1, 0)) %>% 
  ggplot(aes(x= tsne_1, y = tsne_2)) +
  geom_point(aes(colour = pheno_cluster)) #
facet_wrap(~pheno_cluster)

#############################################
library(rstan)

options(mc.cores = parallel::detectCores())
doMC::registerDoMC(cores = 14)
rstan_options(auto_write = TRUE)

model   = stan_model("~/mix_nb_model_Ed_vec_V5.stan")
model   = stan_model("~/Software/cluster_model/mix_nb_model_vec_V10.stan") #mix_nb_model_vec_V6.stan

#vb a function in stan model
clusters = c(25)

runs = c(1:4)

library(doParallel)

foreach(r = runs)%dopar%{
  
  for(number in clusters) {
    
  # number = 25
  # 
  # runs = 1
  
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
  map_est4   = optimizing(model, data = data_cyto, iter = 3000, as_vector = F, init = init_vals0)
  #map_est4   = vb(model, data = data_cyto, algorithm = "fullrank")
  #print(start_time - Sys.time())
  
  direct = paste("Results_by_date/Oct2019/", "full_protein_Ir_init_c_test", number, "_", r, ".rda", sep="") #Desktop/Sep2018/Controls_Oct/Oct_K1/", "K1_mncs_nn_"
  save(map_est4, file = direct)
  
  }
  
  return(start_time - Sys.time())
}

fit_Yule_func = function(cmp_model, data_cyto, num_chains = 4, iter = 1200){
  
  data_for_stan = data_cyto
  
  run_fit = sampling(cmp_model, iter = iter, data = data_for_stan, chains = num_chains, thin = 1, cores = 14)
  
  run_fit
}

run_fit = fit_Yule_func(model, data_cyto)

traceplot(run_fit, inc_warmup =F, pars = "s_rate")
pl1= plot(run_fit, pars = "s_rate")

save(sub_data, file = "Documents/sampled_cell_data_dprior_test_4.rda")

load(file = "Results_by_date/Aug2019/sampled_cell_data_dprior_test_4.rda")

load(file = "Results_by_date/Aug2019/run_dprior_test_4_2.rda") #map_est4 _2

load(file = "Results_by_date/Oct2019/full_protein_Ir_init_c25_4.rda") #map_est4 _2

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

# theta0_c_fit = vals$theta0_s_cp %>% t()

phi_c_fit = vals$phi_s_cp %>% t()

theta1 = vals$theta %>% as.tibble()
colnames(theta1) <- c(1:25)
theta1 %>% mutate(Day = c(1:12)) %>% gather(cluster, prop, -Day) %>% 
  ggplot(aes(x = Day, y = prop))+
  geom_point()+
  facet_wrap(~cluster)
#mu_bg_fit = vals$mu_bg %>% as.vector()
#st_fit = vals$scaling_t
#st_fit = c(1, st_fit)

theta_fit = vals$theta %>% t()
  
  # rep(cluster$vals, length(day_vals)) %>% 
  # matrix(nrow = c)

data_x = gene_data0 %>%  #express_data0 %>%
  filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
  mutate(Day = obstime) %>% 
  #select(Day,  Time) %>% #-cell_cycle, -HBB, -H3, -DNA, -mean_Ir, -mean_Pt ; -cell_id, -cell, -capture, -obstime
  #filter((!Day %in% c("MNCs", "Jurkats"))) %>%
  as.tibble() %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 

data_x1 = data_x %>% select(prot_name) %>% mutate_if(is.numeric, ~ floor(.x)) %>% as.matrix() #Ir191, Ir193, IdU, p_Rb

day_indx = data_x %>%  pull(day_index)

all_cells = predict_cluster(data_x1, theta_fit, mu_c_fit, phi_c_fit, day_indx) # theta0_c_fit
all_cells_prob = predict_cluster_prob(data_x1, theta_fit, mu_c_fit, phi_c_fit, day_indx)

all_cells = x25
x25 = all_cells 
x10 = all_cells

cbind(x25, x10) %>% as.tibble() %>% 
  mutate(c25 = as.factor(V1), c10 = as.factor(V2)) %>% 
  ggplot(aes(x=c25, y =c10))+
  geom_jitter(alpha = 0.5)

# log_prob = calculate_max_prob(data_x1, theta_fit, mu_c_fit, phi_c_fit, day_indx) %>% sum()

loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }


log_like_prob = tibble(no_cluster = c(1:6, 21), log_prob = c(-4862814, -4163492, -4002411, -3961881, -3951784, -3948089, -3924599))

log_like_prob = tibble(no_cluster = c(1:6, 21), log_prob = c(-4862814, -4166917, -4020415, -3981439, -4086365, -3983992, -4314069))

sample_size = log(length(data_x1))

para_size = c * length(prot_list) + c * length(prot_list) + c* length(day_vals)

log_like_prob %>% mutate(BIC_ = -2*log_prob + 
                           sample_size* (no_cluster * length(prot_list) + no_cluster * length(prot_list) + no_cluster* length(day_vals))) %>% 
  ggplot(aes(x = no_cluster, y = BIC_))+
  geom_point()

log_like_prob %>% mutate(ICL_ = log_prob - 
                           0.5* sample_size* (no_cluster * length(prot_list) + no_cluster * length(prot_list) + no_cluster* length(day_vals))) %>% 
  ggplot(aes(x = no_cluster, y = ICL_))+
  geom_point()

#print(table(all_cells)) %>% sum()
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

p2 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(cluster, prob, fill = 0) %>% 
  gather(cluster, prob, -Day) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  ggplot(aes(x = Day, y = prob))+
  geom_point()+
  facet_wrap(~cluster, nrow=2)

p3 = data_x %>% mutate(cluster = all_cells) %>% 
  gather(protein, intensity, CD34:CD41) %>% 
  mutate(Intensity = asinh(intensity)) %>%
  ggplot(aes(x = protein, y = Intensity))+
  geom_violin(size=0.2, scale = "width", trim = FALSE, bw = 0.5)  +
  facet_wrap(~cluster, ncol = 3) 

p4 = signal %>% 
  ggplot(aes(x = prot_name)) +
  geom_point(aes(x = prot_name, y = std, colour = prot_name), size = 3) + #colour = "red"
  # geom_errorbar(aes(x = prot_name, ymax= vals+1.9*std, ymin= vals-1.9*std), width =0.2) +
  theme_classic() +
  theme(text = element_text(size=15), legend.title=element_blank()) +
  facet_wrap(~cluster) +
  scale_y_log10()+
  ylab("Phi") +
  xlab("Protein") 

p4 = signal %>% 
  ggplot(aes(x = prot_name)) +
  geom_point(aes(x = vals, y = std, colour = prot_name), size = 3) + #colour = "red"
  # geom_errorbar(aes(x = prot_name, ymax= vals+1.9*std, ymin= vals-1.9*std), width =0.2) +
  theme_classic() +
  theme(text = element_text(size=15), legend.title=element_blank()) +
  facet_wrap(~cluster) +
  scale_y_log10()+
  scale_x_log10()+
  ylab("Phi") +
  xlab("mu") 

p5 = data_x %>% mutate(cluster = all_cells) %>%
  select(-Time) %>% 
  gather(protein, intensity, mean_Pt) %>% 
  #mutate(Intensity = asinh(intensity)) %>% 
  ggplot()+
  geom_histogram(aes(x = intensity), binwidth = 0.05) +
  scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
  #xlim(0.01, 1000)+
  # geom_vline(xintercept = c(sinh(5), sinh(8)))+
  facet_wrap(~cluster, scales = "free_y")
  facet_grid(cluster~protein, scales = "free")
  
  
############CALCULATE CORRELATIONS
  
cluster = c(1:25)

corr_mean = c()
corr_prop = c()

x_mean = signal %>% select(cluster, prot_name, vals) %>% 
  filter(!prot_name %in% c("IdU", "mean_Ir", "p_Rb"))

x_prop = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(cluster, prob, fill = 0) %>% 
  gather(cluster, prob, -Day) %>% 
  mutate(cluster = as.numeric(cluster)) 

for(i in cluster) {
  
  for(j in cluster){
    
    x1 = x_mean %>% filter(cluster == i) %>% select(vals) %>% t() %>% as.vector()
    x2 = x_mean %>% filter(cluster == j) %>% select(vals) %>% t() %>% as.vector()
    
    test_x = cor.test(x1, x2, method=c("pearson"))
    
    corr_mean_x = test_x$estimate
    
    corr_x = tibble(cluster1 = i, cluster2 = j, corr = corr_mean_x)
    
    corr_mean = rbind(corr_mean, corr_x)
    
  }
  
}

for(i in cluster) {
  
  for(j in cluster){
    
    x1 = x_prop %>% filter(cluster == i) %>% select(prob) %>% t() %>% as.vector()
    x2 = x_prop %>% filter(cluster == j) %>% select(prob) %>% t() %>% as.vector()
    
    test_x = cor.test(x1, x2, method=c("pearson"))
    
    corr_mean_x = test_x$estimate
    
    corr_x = tibble(cluster1 = i, cluster2 = j, corr = corr_mean_x)
    
    corr_prop = rbind(corr_prop, corr_x)
    
  }
  
}

corr_mean %>% mutate(corr_mu = corr, corr_prop = corr_prop$corr) %>%
  mutate(cluster = paste0(cluster1, " ", cluster2)) %>% 
  filter(cluster1 != cluster2) %>% 
  filter(cluster1 <= cluster2) %>% 
  ggplot(aes(x = corr_mu, y = corr_prop))+
  geom_point(aes(colour = as.factor(cluster1))) +
  geom_text(aes(label = cluster), size=4.0) +
  xlim(0.8, 1.0)+
  ylim(0.8, 1.0)

cluster_merge = 
  tibble(cluster = c(1:25),
         merge_cluster = c(1, 2, 3, 4, 5, 1, 7, 8, 9, 1, 11, 12, 13,
                         1, NA, 16, 17, 16, 4, 12, 13, 1, 5, NA, 8))

merge_cluster = c(1, 2, 3, 4, 5, 1, 7, 8, 9, 1, 11, 12, 13,
                  1, 16, 17, 16, 4, 12, 13, 1, 5, 8) %>% unique()
cluster_new = 
  tibble(merge_cluster = merge_cluster,
         state = c(1:length(merge_cluster)))

p1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>%
  left_join(cluster_merge) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>% 
  ggplot(aes(x = prot_name, y= vals)) + 
  geom_point(aes(colour = prot_name)) +
  scale_y_log10()+
  geom_text(aes(label = prot_name), size=3.0)+
  facet_wrap(~state, nrow = 2) +
  theme(legend.position = "none")

p2 = data_x %>% mutate(cluster = as.integer(all_cells)) %>%
  left_join(cluster_merge) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>% 
  group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(state, prob, fill = 0) %>% 
  gather(state, prob, -Day) %>% 
  mutate(state = as.numeric(state)) %>% 
  ggplot(aes(x = Day, y = prob))+
  geom_point()+
  facet_wrap(~state, nrow=2)

gridExtra::grid.arrange(p1, p2, ncol =1)

data1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>%
  left_join(cluster_merge) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) 

data2 = data_x %>% mutate(cluster = as.integer(all_cells)) %>%
  left_join(cluster_merge) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>% 
  group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(state, prob, fill = 0) %>% 
  gather(state, prob, -Day) %>% 
  mutate(state = as.numeric(state))

save(data1, file = "Results_by_date/Dec2019/prot_data.rda")
save(data2, file = "Results_by_date/Dec2019/prop_data.rda")
#############################################

  
  data_x %>% mutate(cluster = all_cells) %>%
    #filter(cluster == 8) %>% #| cluster == 14
    select(-Time) %>% 
    # group_by(obstime) %>%
    # sample_n(5000) %>%
    # ungroup() %>%
    filter(obstime %in% c(6, 8, 16)) %>% 
    select(obstime:capture, mean_Ir, FLI1, IdU, p_Rb, CD38, CD235ab, Pd102:Pd110) %>% 
    # gather(protein, intensity, Pd102:Pd110) %>%
    # mutate(Intensity = intensity, Ir_mean = mean_Ir) %>%
    #mutate(Pd_new = Pd102+Pd106) %>% 
    ggplot()+
    #geom_histogram(aes(x = Intensity), binwidth = 0.05) +
    geom_point(aes(x = CD38, y = Pd102, colour = IdU>sinh(4.5)), alpha = 0.5) +
    #scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
    xlim(0, 200)+
    # geom_vline(xintercept = c(sinh(5), sinh(8)))+
    #facet_wrap(~obstime)
    facet_grid(IdU>sinh(4.5)~obstime) #, scales = "free"
  
  
data_x %>% mutate(cluster = all_cells) %>%
    select(-Time) %>% 
    # group_by(obstime) %>%
    # sample_n(5000) %>%
    # ungroup() %>%
    filter(obstime %in% c(6, 8, 16)) %>% 
    select(obstime:capture, mean_Ir, FLI1, IdU, p_Rb, CD38, CD235ab, Pd102:Pd110) %>% 
    # gather(protein, intensity, Pd102:Pd110) %>%
    # mutate(Intensity = intensity, Ir_mean = mean_Ir) %>%
    #mutate(Pd_new = Pd102+Pd106) %>% 
    ggplot()+
    #geom_histogram(aes(x = Intensity), binwidth = 0.05) +
    geom_point(aes(x = CD38, y = Pd102, colour = IdU>sinh(4.5)), alpha = 0.5) +
    #scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
    xlim(0, 200)+
    # geom_vline(xintercept = c(sinh(5), sinh(8)))+
    #facet_wrap(~obstime)
    facet_grid(IdU>sinh(4.5)~obstime) #, scales = "free"

data_x %>% mutate(cluster = all_cells) %>%
  # filter(Ce140 <= sinh(5.0)) %>% 
  # mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  # filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
  # filter(IdU <= sinh(4)) %>% 
  # group_by(obstime) %>%
  # sample_n(6000) %>%
  # ungroup() %>%
  filter(cluster %in% c(1, 6, 10, 14, 22)) %>% 
  group_by(cluster) %>%
  sample_n(3500) %>%
  ungroup() %>%
  select(-Time) %>% 
  select(obstime:capture, cluster, mean_Ir, mean_Pt, IdU, HBA, CD34, p_Rb, Pd102:Pd110) %>% 
  gather(protein, intensity, Pd102:Pd110) %>% 
  mutate(Intensity = asinh(intensity), Pt_mean = mean_Pt, IdU = IdU) %>% 
  ggplot(aes(x = IdU, y = mean_Ir, colour = asinh(p_Rb)))+
  geom_point(alpha = 0.45) +
  geom_density2d(colour = "red")+
  # facet_grid(cluster~protein)
  facet_wrap(~cluster)

signal %>% 
  filter(cluster %in% c(1, 6, 10)) %>% 
  select(cluster, vals, prot_name) %>% 
  spread(cluster, vals) %>% 
  mutate(ratio1_6 = `1`/`6`, ratio1_10 = `10`/`1`) %>% 
  ggplot(aes(x = prot_name, y = ratio1_6), size = 10) +
  #geom_abline()+
  geom_point(size=5.0) 
  ylim(0, 250)
  facet_grid(~cluster)
  
  
  signal %>% 
    filter(cluster %in% c(1, 14)) %>% 
    select(cluster, vals, prot_name) %>% 
    spread(cluster, vals) %>% 
    #mutate(ratio1_6 = `1`/`6`, ratio1_10 = `10`/`1`) %>% 
    ggplot(aes(x = `1`, y = `14`, colour = prot_name), size = 10) +
    geom_abline()+
    geom_point(size=0.5) +
    geom_text(aes(label = prot_name))+
  ylim(0, 250)
  facet_grid(~cluster)
  
cluster = c(1:25)

pdf(file = paste("Results_by_date/Oct2019/cluster.pdf", sep=""), width = 15, height = 10) 

for(i in cluster){
  ref_data = signal %>% select(cluster, vals, prot_name) %>% 
  filter(cluster == i) %>% 
  mutate(ref_vals = vals) %>% 
  select(-vals)

  p1 = signal %>% left_join(ref_data, by = c("prot_name")) %>%
  ggplot(aes(x = vals, y = ref_vals, colour = prot_name), size = 1.5) +
  geom_abline()+
  geom_point(size=0.2) +
  geom_text(aes(label = prot_name))+
  theme(legend.position = "none")+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~cluster.x)+
  ggtitle(i)
  
  plot(p1)
  #pdf(file = paste("Results_by_date/Oct2019/", i, ".pdf", sep="")) 
  
}

  dev.off() 
  
  
data_x %>% mutate(cluster = all_cells) %>%
  filter(cluster %in% c(15)) %>% 
  select(-Time) %>%

asinh(20)
sinh(4)
#####*  
pro_data = gene_data0 %>%
    filter(Ce140 <= sinh(5.0)) %>% 
    mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
    filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

d1 = pro_data %>% filter(obstime == 0) %>% 
    select(-Time) %>% 
    select(obstime:capture, Pd102:Pd105) %>%
    mutate(mean_Pd0 = mean_irid(Pd102, Pd104)) %>% 
    mutate(mean_Pd = mean_irid(mean_Pd0, Pd105)) %>% 
    select(obstime:capture, mean_Pd)

d2 = pro_data %>% filter(obstime == 2) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd104, Pd106) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd104)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd106)) %>% 
  select(obstime:capture, mean_Pd)

d3 = pro_data %>% filter(obstime == 4) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd104, Pd108) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd104)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd108)) %>% 
  select(obstime:capture, mean_Pd)

d4 = pro_data %>% filter(obstime == 6) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd104, Pd110) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd104)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd110)) %>% 
  select(obstime:capture, mean_Pd)

d5 = pro_data %>% filter(obstime == 8) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd105, Pd106) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd105)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd106)) %>% 
  select(obstime:capture, mean_Pd)
  
d6 = pro_data %>% filter(obstime == 10) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd105, Pd108) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd105)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd108)) %>% 
  select(obstime:capture, mean_Pd)

d7 = pro_data %>% filter(obstime == 11) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd105, Pd110) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd105)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd110)) %>% 
  select(obstime:capture, mean_Pd)
  
d8 = pro_data %>% filter(obstime == 12) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd105, Pd108) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd105)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd108)) %>% 
  select(obstime:capture, mean_Pd)

d9 = pro_data %>% filter(obstime == 14) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd106, Pd110) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd106)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd110)) %>% 
  select(obstime:capture, mean_Pd)

d10 = pro_data %>% filter(obstime == 16) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd102, Pd108, Pd110) %>%
  mutate(mean_Pd0 = mean_irid(Pd102, Pd108)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd110)) %>% 
  select(obstime:capture, mean_Pd)

d11 = pro_data %>% filter(obstime == 18) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd104, Pd106, Pd110) %>%
  mutate(mean_Pd0 = mean_irid(Pd104, Pd106)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd110)) %>% 
  select(obstime:capture, mean_Pd)

d12 = pro_data %>% filter(obstime == 20) %>% 
  select(-Time) %>% 
  select(obstime:capture, Pd104, Pd108, Pd110) %>%
  mutate(mean_Pd0 = mean_irid(Pd104, Pd108)) %>% 
  mutate(mean_Pd = mean_irid(mean_Pd0, Pd110)) %>% 
  select(obstime:capture, mean_Pd)

d1 %>% rbind(d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12) %>% 
  ggplot(aes(x = obstime, y = asinh(mean_Pd), fill = obstime))+
  geom_violin(size=0.2, scale = "width", trim = FALSE, bw = 0.5)

gene_data0 %>%
  filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>%
  select(-Time) %>% 
  select(obstime:capture, ATRX, CD33, mean_Pt, mean_Ir, FLI1) %>% 
  gather(protein, intensity, ATRX:FLI1) %>% 
  mutate(Intensity = asinh(intensity)) %>% 
  ggplot(aes(x = obstime, y = Intensity, fill = obstime))+
  geom_violin(size=0.2, scale = "width", trim = FALSE, bw = 0.5)+
  facet_wrap(~protein, scales = "free")

#####FIND CLUSTER 1
gene_data0 %>%
  filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
  select(-Time) %>% 
  filter(IdU<=sinh(3.5)) %>% 
  #select(obstime:capture, mean_Ir, mean_Pt, IdU, CD34, HBA, Pd102:Pd110) %>% 
  gather(protein, intensity, Pd102:Pd110) %>%
  filter(intensity <400, mean_Ir >= 90, mean_Ir <= 200) %>% #, intensity <400) %>% 
  select(obstime:capture, ATRX:mean_Ir) %>% 
  gather(protein, intensity, ATRX:mean_Ir) %>%
  mutate(Intensity = asinh(intensity)) %>% 
  ggplot(aes(x = obstime, y = Intensity, fill = obstime))+
  geom_violin(size=0.2, scale = "width", trim = FALSE, bw = 0.5)+
  facet_wrap(~protein, scales = "free")

gene_data0 %>%
    filter(Ce140 <= sinh(5.0)) %>% 
    mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
    filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
    filter(IdU <= sinh(4)) %>% 
    group_by(obstime) %>% 
    sample_n(2000) %>% 
    ungroup() %>%
    filter(obstime %in% c(0, 4, 6, 14, 16)) %>% 
    select(-Time) %>% 
    select(obstime:capture, mean_Ir, mean_Pt, IdU, CD34, HBA, Pd102:Pd110) %>% 
    gather(protein, intensity, Pd102:Pd110) %>% 
    mutate(Intensity = asinh(intensity), Pt_mean = mean_Pt, IdU = asinh(IdU)) %>% 
    ggplot()+
    #geom_histogram(aes(x = Intensity), binwidth = 0.05) +
    geom_point(aes(x = intensity, y = mean_Ir, colour = IdU), alpha = 0.45) +
    #geom_contour(aes(x = Intensity, y = Ir_mean, z = Ir_mean))+
    #scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
    #ylim(0, 50)+
    # geom_vline(xintercept = c(sinh(5), sinh(8)))+
    # facet_wrap(~obstime)
    facet_grid(obstime~protein) 
  
  gene_data0 %>%
    filter(Ce140 <= sinh(5.0)) %>% 
    mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
    filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
    group_by(obstime) %>% 
    sample_n(2000) %>% 
    ungroup() %>%
    select(-Time) %>% 
    select(obstime:capture, IdU, mean_Ir, Pd102:Pd110) %>%
    gather(protein, intensity, Pd102:Pd110) %>%
    mutate(IdU_l = ifelse(IdU>=sinh(4.5), T, F)) %>% 
    mutate(new_Ir = ifelse(intensity <= 27.28992, NA, mean_Ir/(intensity))) %>% 
    # filter(obstime == 4 | obstime == 10 | obstime == 20) %>% 
    filter(new_Ir<=5) %>% 
    ggplot()+
    geom_histogram(aes(x = asinh(mean_Ir)), bins = 100) + #, binwidth = 0.05
    # geom_point(aes(x = asinh(Pd106), y = asinh(FLI1)), alpha = 0.45) +
    #geom_contour(aes(x = Intensity, y = Ir_mean, z = Ir_mean))+
    #scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
    #xlim(0.01, 650)+
    # geom_vline(xintercept = c(sinh(5), sinh(8)))+
    # facet_wrap(~obstime, scales = "free")
    facet_grid(IdU_l~obstime, scales = "free") #, scales = "free"
  
  gene_data0 %>%
    filter(Ce140 <= sinh(5.0)) %>% 
    mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
    filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
    group_by(obstime) %>% 
    sample_n(2000) %>% 
    ungroup() %>%
    select(-Time) %>% 
    select(obstime:capture, mean_Ir, Pd102:Pd110) %>%
    gather(protein, intensity, Pd102:Pd110) %>%
    mutate(new_Ir = ifelse(intensity <= 27.28992, NA, mean_Ir/(intensity))) %>% 
    # filter(obstime == 4 | obstime == 10 | obstime == 20) %>% 
    filter(new_Ir>=3, new_Ir<=5) %>% 
    ggplot()+
    geom_histogram(aes(x = new_Ir), bins = 100) +
    facet_grid(obstime~protein, scales = "free") 
    
  gene_data0 %>%
      filter(Ce140 <= sinh(5.0)) %>% 
      mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
      group_by(obstime) %>% 
      sample_n(2000) %>% 
      ungroup() %>%
      select(-Time) %>% 
      #select(obstime:capture, mean_Ir, IdU, Pd102:Pd110) %>% 
      #filter(Pd104 >= 200 & Pd106 >= 200) %>% 
      #gather(protein, intensity, Pd102:Pd110) %>% 
      #mutate(Intensity = asinh(intensity), Ir_mean = asinh(mean_Ir), IdU = asinh(IdU)) %>% 
      ggplot()+
      #geom_histogram(aes(x = Intensity), binwidth = 0.05) +
      geom_point(aes(x = asinh(Pd106), y = asinh(FLI1)), alpha = 0.45) +
      #geom_contour(aes(x = Intensity, y = Ir_mean, z = Ir_mean))+
      #scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
      #ylim(0.01, 650)+
      # geom_vline(xintercept = c(sinh(5), sinh(8)))+
      facet_wrap(~obstime, scales = "free")
    facet_grid(obstime~protein) #, scales = "free"

  data_low = gene_data0 %>%
    filter(Ce140 <= sinh(5.0)) %>% 
    mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
    filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
    group_by(obstime) %>% 
    sample_n(1000) %>% 
    ungroup() %>%
    select(-Time) %>% 
    #select(obstime:capture, mean_Ir, IdU, Pd102:Pd110) %>% 
    filter(Pd104 < 200 & Pd106 < 200) %>% 
    gather(protein, intensity, ATRX:mean_Ir) %>% 
    filter(obstime==2) %>% 
    # mutate(Intensity = asinh(intensity)) %>% 
    group_by(protein) %>% 
    summarise(mu_in = mean(intensity)) %>% 
    mutate(cat = "Low")
  
  data_high = gene_data0 %>%
    filter(Ce140 <= sinh(5.0)) %>% 
    mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
    filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
    group_by(obstime) %>% 
    sample_n(1000) %>% 
    ungroup() %>%
    select(-Time) %>% 
    #select(obstime:capture, mean_Ir, IdU, Pd102:Pd110) %>% 
    filter(Pd104 >= 200 & Pd106 >= 200) %>% 
    gather(protein, intensity, ATRX:mean_Ir) %>% 
    filter(obstime==2) %>% 
    # mutate(Intensity = asinh(intensity)) %>% 
    group_by(protein) %>% 
    summarise(mu_in = mean(intensity)) %>% 
    mutate(cat = "High")
  
  rbind(data_low, data_high) %>% 
    spread(cat, mu_in) %>% 
    filter(!protein %in% c("Pd102", "Pd104", "Pd105", "Pd106", "Pd108", "Pd110")) %>% 
    ggplot()+
    geom_hline(yintercept = 1)+
    geom_point(aes(x = protein, y = High/Low)) 
    
  gene_data0 %>%
      filter(Ce140 <= sinh(5.0)) %>% 
      mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
      filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
      group_by(obstime) %>% 
      sample_n(1000) %>% 
      ungroup() %>%
      select(-Time) %>% 
      #select(obstime:capture, mean_Ir, IdU, Pd102:Pd110) %>% 
      filter(Pd104 < 200 & Pd106 < 200) %>% 
      gather(protein, intensity, ATRX:mean_Ir) %>% 
      filter(obstime==2) %>% 
      # mutate(Intensity = asinh(intensity)) %>% 
      group_by(protein) %>% 
      summarise(mu_in = mean(intensity)) %>%
      ggplot()+
      # geom_histogram(aes(x = Intensity), binwidth = 0.05) +
      geom_point(aes(x = protein, y = mu_in)) #+
    #geom_contour(aes(x = Intensity, y = Ir_mean, z = Ir_mean))+
    #scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
    #ylim(0.01, 650)+
    xlim(0.01, 8)+
    # geom_vline(xintercept = c(sinh(5), sinh(8)))+
     facet_wrap(~protein, scales = "free")
    facet_grid(obstime~protein) #, scales = "free"
  
png("Documents/file.png",width=6,height=4) 
dev.off()
  
  data_x %>% mutate(cluster = all_cells) %>%
    #filter(cluster == 8) %>% #| cluster == 14
    select(-Time) %>% 
    select(obstime:capture, cluster, Pd102:Pd110) %>% 
    gather(protein, intensity, Pd102:Pd110) %>% 
    mutate(Intensity = asinh(intensity)) %>% 
    ggplot()+
    geom_histogram(aes(x = Intensity), binwidth = 0.05) +
    #scale_x_continuous(trans = 'asinh', limits = c(0.01, 1000))+
    #xlim(0.01, 1000)+
    # geom_vline(xintercept = c(sinh(5), sinh(8)))+
    #facet_wrap(~cluster, scales = "free_y")
    facet_grid(obstime~protein, scales = "free")

all_cells0 <- all_cells[,1]

p6 = data_x %>% mutate(cluster = all_cells0) %>%
  mutate(state = cluster) %>% 
  ggplot()+
  geom_point(aes(x = asinh(CD235ab), y = asinh(HBA), colour = asinh(CD34)), alpha = 0.2)+
  facet_wrap(~cluster)

data_x %>% mutate(cluster = all_cells0) %>%
  left_join(cluster_order, by = c("cluster" = "old_order")) %>% 
  select(Day, IdU,new_order) %>%
  # mutate(Intensity = asinh(intensity)) %>%
  #filter(new_order == 2) %>% 
  # gather(protein, intensity, CD34:CD41) %>% 
  # mutate(Intensity = asinh(intensity)) %>%
  # ggplot(aes(x = protein, y = Intensity, fill = protein))+
  # geom_violin(size=0.2, scale = "width", trim = FALSE, bw = 0.5)  +
  # facet_wrap(~Day, ncol = 3) 
  ggplot(aes(x = IdU))+
  geom_histogram(bins = 40)+
  #geom_point(alpha = 0.15)+
  #scale_y_log10()+
  scale_x_continuous(trans = 'asinh')+
  #scale_y_continuous(trans = 'asinh')+
  facet_grid(new_order~Day) #, scales = "free_y"

data_x %>% mutate(cluster = all_cells0) %>%
  #left_join(cluster_order, by = c("cluster" = "old_order")) %>% 
  select(Day, IdU,new_order) %>%
  mutate(Cell_cycle = ifelse(IdU > sinh(3), "Cycle", "Non-cycle"))
  # mutate(Intensity = asinh(intensity)) %>%
  #filter(new_order == 2) %>% 
  # gather(protein, intensity, CD34:CD41) %>% 
  # mutate(Intensity = asinh(intensity)) %>%
  # ggplot(aes(x = protein, y = Intensity, fill = protein))+
  # geom_violin(size=0.2, scale = "width", trim = FALSE, bw = 0.5)  +
  # facet_wrap(~Day, ncol = 3) 
  ggplot(aes(x = IdU))+
  geom_histogram(bins = 40)+
  #geom_point(alpha = 0.15)+
  #scale_y_log10()+
  scale_x_continuous(trans = 'asinh')+
  #scale_y_continuous(trans = 'asinh')+
  facet_grid(new_order~Day) #, scales = "free_y"
  
  ########Calculate S phase##################

  S_Phase_ratio = data_x %>% mutate(cluster = all_cells0) %>%
    left_join(cluster_order, by = c("cluster" = "old_order")) %>% 
    select(Day, IdU, new_order) %>%
    mutate(S_Phase = ifelse(IdU <= sinh(3.466711), "Other", "S_Phase")) %>% 
    group_by(new_order) %>% 
    mutate(total = n()) %>% 
    ungroup() %>% 
    group_by(new_order, S_Phase) %>% 
    mutate(count = n()/total) %>% ungroup() %>% 
    select(new_order, S_Phase, count) %>% 
    unique() %>% 
    filter(S_Phase == "S_Phase") %>% 
    arrange(new_order) %>% 
    select(-S_Phase) %>%
    pull(count) %>% 
    as.vector()
    
    ggplot()+
    geom_bar(aes(x = Day, y = (..count..)/sum(..count..), fill = S_Phase)) +
    #scale_y_log10()+
    facet_grid(new_order~Day)
  
  data_x %>% mutate(cluster = all_cells0) %>%
    left_join(cluster_order, by = c("cluster" = "old_order")) %>% 
    select(Day, IdU, CD235ab, HBA, new_order) %>%
    mutate(Cell_cycle = ifelse(IdU <= sinh(3.5), "Others", "S Phase")) %>% 
    group_by(Day) %>% 
    mutate(tot_count = n()) %>% ungroup() %>% 
    group_by(Day, Cell_cycle) %>% 
    mutate(count = n()/tot_count) %>% 
    ungroup() %>% 
    select(Day, Cell_cycle, count, new_order) %>% 
    unique() %>% 
    ggplot()+
    geom_bar(aes(x = Cell_cycle, y =count, fill = Cell_cycle), stat="identity") +
    #scale_y_log10()+
    facet_grid(new_order~Day)


gridExtra::grid.arrange(gridExtra::grid.arrange(p1, p2, ncol =1), p5, heights = c(2, 0.5), ncol=1)

vals$mu_s_cp %>% as.tibble() %>% 
  mutate(cluster = 1:c) %>% 
  gather(protein, vals, -cluster) %>%
  mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
  left_join(prot_list, by = c("protein" = "prot_index")) %>% 
  left_join(phi) %>% 
  rnbinom( 100, 1, 1)

#################simulation data#################
pred_data = data_x %>% mutate(cluster = all_cells) %>%
  select(-Time) %>% 
  gather(protein, intensity, CD34:CD235ab)

pred_model = vals$mu_s_cp %>% as.tibble() %>% 
  mutate(cluster = 1:c) %>% 
  gather(protein, vals, -cluster) %>%
  mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
  left_join(prot_list, by = c("protein" = "prot_index")) %>% 
  left_join(phi) %>% 
  select(-protein)

pred_cluster = pred_data %>% select(cluster) %>% 
  unique() %>% arrange(cluster) %>% pull(cluster)

pred_prot =  pred_data %>% select(protein) %>% 
  unique() %>% arrange(protein) %>% pull(protein)

# pp_out = tibble(x = data_x) %>% 
#   ggplot(aes(x)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y), col = "red", alpha= 0.9) +
#   ggtitle(protein) + ylab("Proportion of counts") #+  xlim(-1, 50)
# list(pp = pp_out, fit_out = fit_out)

tbl_pred = c()

for(i in pred_cluster){
  
  for(j in pred_prot) {
    
    # i = 2
    # j = "CD34"
    print(i)
    
    data_x_select = pred_data %>% as.tibble() %>% dplyr::filter(cluster == i) %>% dplyr::filter(protein == j) %>% pull(intensity)
    
    mu = pred_model %>% as.tibble() %>% dplyr::filter(cluster == i) %>% dplyr::filter(prot_name  == j) %>% pull(vals)
    
    phi = pred_model %>% as.tibble() %>% dplyr::filter(cluster == i) %>% dplyr::filter(prot_name  == j) %>% pull(std)
    
    tbl_pred = tbl_pred %>% rbind(tibble(cluster = i, protein = j, 
                      x = 0:max(data_x_select), y = dnbinom(x, mu = mu, size = phi)))
    
  }
}

plot_data = data_x %>% mutate(cluster = all_cells) %>%
  select(-Time) %>% 
  gather(protein, intensity, CD34:CD235ab) %>% 
  group_by(protein, cluster, intensity) %>% 
  # summarise(count = n()) %>% 
  summarise(count = n()) %>%
  ungroup() %>% 
  mutate(count = as.numeric(count)) %>%
  group_by(protein, cluster) %>% 
  mutate(perc = sum(count)) %>% 
  ungroup() %>% 
  mutate(perc0 = count/perc) 

plot_data %>% 
  as.tibble() %>% 
  ggplot(aes(x = intensity, y=perc0))+
  geom_col(width = 0.2) +
  #geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +
  geom_point(data = tbl_pred, aes(x = x, y =y), col = "red", alpha= 0.6, size = 1.0)+
  scale_x_continuous(trans = 'asinh')+
  facet_grid(cluster~protein, scales = "free")

total_n = 7500

# pred_model %>% View()

cluster1 = tibble(cluster = 1,
                  CD34 = rnbinom(total_n, mu = 59.7303771, size = 0.9118761), 
                  CD235ab = rnbinom(total_n, mu = 5.6809071, size = 1.5420040),
                  HBA = rnbinom(total_n, mu = 0.2399894, size = 0.1631683))

cluster2 = tibble(cluster = 2,
                  CD34 = rnbinom(total_n, mu = 0.8153499, size = 0.1713387), 
                  CD235ab = rnbinom(total_n, mu = 5.9681089, size = 2.3663730),
                  HBA = rnbinom(total_n, mu = 0.1514670, size = 0.1651301))

cluster3 = tibble(cluster = 3,
                  CD34 = rnbinom(total_n, mu = 0.3121052, size = 0.1688198), 
                  CD235ab = rnbinom(total_n, mu = 228.9127665, size = 4.2789123),
                  HBA = rnbinom(total_n, mu = 2.2062987, size = 0.2111989))

cluster4 = tibble(cluster = 4,
                  CD34 = rnbinom(total_n, mu = 0.4460181, size = 0.1538042), 
                  CD235ab = rnbinom(total_n, mu = 223.9129869, size = 6.8593690),
                  HBA = rnbinom(total_n, mu = 1395.4081914, size = 1.4528272))

cluster1 %>% rbind(cluster2, cluster3, cluster4) %>% 
  gather(protein, intensity, CD34:HBA) %>% 
  ggplot(aes(x = intensity))+
 geom_histogram(bins = 50)+
  scale_x_continuous(trans = 'asinh')+
  facet_grid(cluster~protein, scales = "free")
  
simu_vals = cluster1 %>% rbind(cluster2, cluster3, cluster4)

gated_data = simu_vals

prot_name = gated_data %>% select(-cluster) %>% colnames() %>% as.character()

prot_list = data.frame(prot_name = prot_name, prot_index = 1:length(prot_name)) %>% as.tibble() %>% 
  mutate(prot_name = as.character(prot_name))

#summarise the gated_data
sub_data = gated_data 

x_data   = sub_data %>% select(-cluster) %>% as.matrix()
vec_d    = sub_data %>% mutate(new_indx = as.numeric(factor(as.character(cluster), levels = unique(as.numeric(cluster))))) %>% pull(new_indx)
indx_tp  = sub_data %>% group_by(cluster) %>% summarise(num_cells = n()) %>% mutate(num_cells  = cumsum(num_cells)) %>% 
  mutate(from = c(0, num_cells[-n()]) + 1) %>% select(from, num_cells) %>% as.matrix()

x_indx = x_data %>% as_tibble() %>% mutate_all(~as.numeric(as.factor(.x)))
x_unique_aux = x_data %>% as_tibble() %>% map2(1:ncol(.), ~ as.factor(.x) %>% levels() %>% 
                                                 tibble(unique_count = ., prot =.y)) %>% bind_rows()
unique_prot_x = x_unique_aux %>% mutate(unique_count = as.numeric(unique_count)) %>% pull(unique_count)
indx_prot  = x_unique_aux %>% group_by(prot) %>% summarise(num_vals = n()) %>% mutate(num_vals  = cumsum(num_vals)) %>% 
  mutate(from = c(0, num_vals[-n()]) + 1) %>% select(from, num_vals) %>% as.matrix()

data_cyto = list(c = 4, k = ncol(x_data), N = nrow(x_data), x = x_data, d = max(vec_d), 
                 day_index = vec_d, indx_tp = indx_tp,
                 num_unique_counts = length(unique_prot_x),
                 unique_prot_x = unique_prot_x,
                 indx_prot = indx_prot,
                 x_indx = x_indx)

vec_d %>% unique()

library(rstan)

options(mc.cores = parallel::detectCores())
doMC::registerDoMC(cores = 14)
rstan_options(auto_write = TRUE)

model   = stan_model("~/mix_nb_model_Ed_vec_V5.stan")
model   = stan_model("~/Software/cluster_model/mix_nb_model_vec_V8.stan") #mix_nb_model_vec_V6.stan

#########################################
model   = stan_model("~/Software/cluster_model/mix_nb_model_vec_V10.stan") #mix_nb_model_vec_V6.stan

#vb a function in stan model
runs = c(1:10)

library(doParallel)

foreach(number = runs)%dopar%{
  
  print(number)
  
  data_cyto = list(c = number, k = ncol(x_data), N = nrow(x_data), x = x_data, d = max(vec_d), 
                   day_index = vec_d, indx_tp = indx_tp,
                   num_unique_counts = length(unique_prot_x),
                   unique_prot_x = unique_prot_x,
                   indx_prot = indx_prot,
                   x_indx = x_indx)
  
  start_time <- Sys.time()
  map_est4   = optimizing(model, data = data_cyto, iter = 5000, as_vector = F)
  #map_est4   = vb(model, data = data_cyto, algorithm = "fullrank")
  #print(start_time - Sys.time())
  
  direct = paste("Results_by_date/Sep2019/", "simulate_density_prob_", number, ".rda", sep="") #Desktop/Sep2018/Controls_Oct/Oct_K1/", "K1_mncs_nn_"
  save(map_est4, file = direct)
  
  return(start_time - Sys.time())
}
#######################################
#vb a function in stan model
runs = c(1:4)

library(doParallel)

foreach(number = runs)%dopar%{
  start_time <- Sys.time()
  map_est4   = optimizing(model, data = data_cyto, iter = 5000, as_vector = F)
  #map_est4   = vb(model, data = data_cyto, algorithm = "fullrank")
  #print(start_time - Sys.time())
  
  direct = paste("Results_by_date/Sep2019/", "run_dprior_test_simu_ng_00", number, ".rda", sep="") #Desktop/Sep2018/Controls_Oct/Oct_K1/", "K1_mncs_nn_"
  save(map_est4, file = direct)
  
  return(start_time - Sys.time())
}

save(sub_data, file = "Results_by_date/Aug2019/sampled_cell_data_dprior_test_4.rda")

load(file = "Results_by_date/Aug2019/sampled_cell_data_dprior_test_4.rda")

load(file = "Results_by_date/Aug2019/run_dprior_test_4_2.rda") #map_est4 _2

load(file = "Results_by_date/Sep2019/simulate_density_prob_4.rda") #map_est4 _2

run_fit = map_est4

day_vals = sub_data$cluster %>% unique()
day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)

c = 4

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

day_vals = gated_data$cluster %>% unique()
day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)

mu_c_fit = vals$mu_s_cp %>% t()

theta0_c_fit = vals$theta0_s_cp %>% t()

phi_c_fit = vals$phi_s_cp %>% t()

#mu_bg_fit = vals$mu_bg %>% as.vector()
#st_fit = vals$scaling_t
#st_fit = c(1, st_fit)

theta_fit = vals$theta %>% t()

# rep(cluster$vals, length(day_vals)) %>% 
# matrix(nrow = c)

data_x = gated_data %>%  #express_data0 %>%
  mutate(Day = cluster) %>% 
  select(Day, CD34, HBA, CD235ab, mean_Ir) %>% #-cell_cycle, -HBB, -H3, -DNA, -mean_Ir, -mean_Pt ; -cell_id, -cell, -capture, -obstime
  mutate_if(is.numeric, ~ floor(.x)) %>% 
  #filter((!Day %in% c("MNCs", "Jurkats"))) %>%
  as.tibble() %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 

data_x1 = data_x %>% select(-Day, -day_index) %>% as.matrix() #Ir191, Ir193, IdU, p_Rb

day_indx = data_x %>%  pull(day_index)

all_cells = predict_cluster(data_x1, theta_fit, mu_c_fit, phi_c_fit, day_indx) # 

#print(table(all_cells)) %>% sum()
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

p2 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(cluster, prob, fill = 0) %>% 
  gather(cluster, prob, -Day) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  ggplot(aes(x = Day, y = prob))+
  geom_point()+
  facet_wrap(~cluster, nrow=2)

gridExtra::grid.arrange(p1, p2, ncol =1)

p5 = data_x %>% mutate(cluster = all_cells) %>%
  gather(protein, intensity, CD34:mean_Ir) %>% 
  #mutate(Intensity = asinh(intensity)) %>% 
  ggplot()+
  geom_histogram(aes(x = intensity), binwidth = 0.5) +
  #scale_x_continuous(trans = 'asinh')+
  #xlim(0.01, 9)+
  facet_grid(cluster~protein, scales = "free")

#################simulation data#################
pred_data = data_x %>% mutate(cluster = all_cells) %>%
  #select(-Time) %>% 
  gather(protein, intensity, CD34:CD235ab)

pred_model = vals$mu_s_cp %>% as.tibble() %>% 
  mutate(cluster = 1:c) %>% 
  gather(protein, vals, -cluster) %>%
  mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
  left_join(prot_list, by = c("protein" = "prot_index")) %>% 
  left_join(phi) %>% 
  select(-protein)

pred_cluster = pred_data %>% select(cluster) %>% 
  unique() %>% arrange(cluster) %>% pull(cluster)

pred_prot =  pred_data %>% select(protein) %>% 
  unique() %>% arrange(protein) %>% pull(protein)

# pp_out = tibble(x = data_x) %>% 
#   ggplot(aes(x)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y), col = "red", alpha= 0.9) +
#   ggtitle(protein) + ylab("Proportion of counts") #+  xlim(-1, 50)
# list(pp = pp_out, fit_out = fit_out)

tbl_pred = c()

for(i in pred_cluster){
  
  for(j in pred_prot) {
    
    # i = 2
    # j = "CD34"
    print(i)
    
    data_x_select = pred_data %>% as.tibble() %>% dplyr::filter(cluster == i) %>% dplyr::filter(protein == j) %>% pull(intensity)
    
    mu = pred_model %>% as.tibble() %>% dplyr::filter(cluster == i) %>% dplyr::filter(prot_name  == j) %>% pull(vals)
    
    phi = pred_model %>% as.tibble() %>% dplyr::filter(cluster == i) %>% dplyr::filter(prot_name  == j) %>% pull(std)
    
    tbl_pred = tbl_pred %>% rbind(tibble(cluster = i, protein = j, 
                                         x = 0:max(data_x_select), y = dnbinom(x, mu = mu, size = phi)))
    
  }
}

plot_data = data_x %>% mutate(cluster = all_cells) %>%
  #select(-Time) %>% 
  gather(protein, intensity, CD34:CD235ab) %>% 
  group_by(protein, cluster, intensity) %>% 
  # summarise(count = n()) %>% 
  summarise(count = n()) %>%
  ungroup() %>% 
  mutate(count = as.numeric(count)) %>%
  group_by(protein, cluster) %>% 
  mutate(perc = sum(count)) %>% 
  ungroup() %>% 
  mutate(perc0 = count/perc) 

plot_data %>% 
  as.tibble() %>% 
  ggplot(aes(x = intensity, y=perc0))+
  geom_col(width = 0.2) +
  #geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +
  geom_point(data = tbl_pred, aes(x = x, y =y), col = "red", alpha= 0.6, size = 1.0)+
  scale_x_continuous(trans = 'asinh')+
  facet_grid(cluster~protein, scales = "free")
#############################################
library(rstan)
library(tidyverse)
############dirichlet fit############  
cmp_model = stan_model(file = "Results_by_date/Aug2019/State_transitions_w_NetGrowth_test_v3.stan") #Results_by_date/Aug2019/State_transitions_w_NetGrowth.stan

cluster_order = tibble(new_order = c(2,1,3,4), old_order = c(1:4))

data_t = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
  filter(Day != 0, Day != 11) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(cluster, prob, fill = 0) %>% 
  gather(cluster, prob, -Day) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  left_join(cluster_order, by = c("cluster" = "old_order")) %>% 
  select(-cluster) %>% 
  spread(new_order, prob) %>% 
  select(-Day) %>% 
  as.matrix()

#growth_total = c(1, 2.1, 4.8, 4.2, 3.1, 5.7, 1.1, 1.7, 1.5, 0.7, 0.7, 0.8) * 7 * 1e6 #1.1,

S_Phase_ratio = data_x %>% mutate(cluster = all_cells0) %>%
  left_join(cluster_order, by = c("cluster" = "old_order")) %>% 
  select(Day, IdU, new_order) %>%
  mutate(S_Phase = ifelse(IdU <= sinh(3.466711), "Other", "S_Phase")) %>% 
  group_by(new_order) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  group_by(new_order, S_Phase) %>% 
  mutate(count = n()/total) %>% ungroup() %>% 
  select(new_order, S_Phase, count) %>% 
  unique() %>% 
  filter(S_Phase == "S_Phase") %>% 
  arrange(new_order) %>% 
  select(-S_Phase) %>%
  pull(count) %>% 
  as.vector()


growth_total = c(14.7, 70.42, 295.76, 921.1, 5205.06,  9553.36, 14610.96, 11044.95, 8242.5, 6594.) * 1e6 #1.1,7., 5619.63,

# growth_total %>% as.vector() %>% as.integer()
data_y = growth_total %>% ceiling()
typeof(data_y)
knn_list = tibble(state.x = c(2, 3, 4),
                  state.y = c(1, 2, 3))

n_state = ncol(data_t)
n_time = nrow(data_t)
n_rate = nrow(knn_list)

time_points = data_x %>% select(Day) %>% filter(Day !=0, Day != 11) %>% unique() %>% pull(Day) %>% as.integer()

gen_init_vals = function(n_rate, n_state, n_time){
  list(#lambda1 = runif(n = n_rate, 1, 2),
       #lambda_net = runif(n = n_rate*n_state, -2, 2),
       Sigma_vec1 = runif(n = n_state, 1000, 3000),
       Sigma_vec2 = runif(n = n_state, 1000, 3000),
       Sigma_vec3 = runif(n = n_state, 1000, 3000),
       lambda_trans1 = runif(n = 1, 0, 2),
       lambda_trans2 = runif(n = 1, 0, 1),
       Sigma_tot = runif(n = 1, 0, 1))
} # true_y = runif(n = n_time, 1e5, 5e5)

fit_Yule_func = function(cmp_model, data_t, data_y, n_state, n_time, time_points, n_rate, knn_list, S_Phase_ratio, num_chains = 4, iter = 1000){
  
  n_state = ncol(data_t)
  n_time = nrow(data_t)
  n_rate = nrow(knn_list)
  
  data_for_stan = list(n_state = n_state, n_time = n_time, x_t_in = data_t, y_t_in = data_y, #y_t_sd_in = data_y_sd, 
                       time_points = time_points, n_rate_in=n_rate, knn_weight = knn_list,
                       s_phase_ratio = S_Phase_ratio) 
  
  init_vals = rerun(num_chains, gen_init_vals(n_rate, n_state, n_time))
  
  run_fit = sampling(cmp_model, iter = iter, data = data_for_stan, chains = num_chains, init = init_vals, thin = 1, cores = 14) # warmup = 500, init = init_vals
  
  run_fit
}
#, init = init_vals
# fit_Yule_func = function(cmp_model, data_t, n_state, n_time, time_points, n_rate, knn_list){ , init = init_vals
#   
#   num_chains = 4
#   
#   data_for_stan = list(n_state = n_state, n_time = n_time, x_t_in = data_t, time_points = time_points, n_rate_in=n_rate, knn_weight = knn_list) 
#     
#   
#   run_fit = sampling(cmp_model, iter = 1000, warmup = 500, data = data_for_stan, chains = num_chains, thin = 1, cores = 14)
#   
#   run_fit
# }

run_fit = fit_Yule_func(cmp_model, data_t, data_y, n_state, n_time, time_points, n_rate, knn_list, S_Phase_ratio)

save(run_fit, file = "Results_by_date/Aug2019/3_protein_rates_trans.rda")
load(file = "Results_by_date/Aug2019/3_protein_rates_trans.rda")

traceplot(run_fit, inc_warmup =F, pars = "s_rate")
pl1= plot(run_fit, pars = "s_rate")
pl2= plot(run_fit,  pars = "lambda_trans1")
pl3= plot(run_fit,  pars = "lambda_trans2")
pl4= plot(run_fit,  pars = "lambda_death")

p1 = plot(run_fit, pars = "lambda1")
p2 = plot(run_fit, pars = "lambda2")
p3 = plot(run_fit, pars = "lambda3")

pp1 = gridExtra::grid.arrange(pl1, pl2, pl3, pl4, nrow = 4)

tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
tbl_summary %>% View()
alpha1 = rstan::summary(run_fit, par = "alpha1")
alpha2 = rstan::summary(run_fit, par = "alpha2")
alpha3 = rstan::summary(run_fit, par = "alpha3")

##########################################################
#########################################################################
new_time_points = seq(2, 20, 0.1)
new_time_length = length(new_time_points)
new_time_index = tibble(time = new_time_points, time_points = 1:new_time_length)

x_t_1 = rstan::summary(run_fit, par = "x_t_1") 

x_t_1_new = x_t_1$summary  

predict = x_t_1_new %>% as.tibble(rownames = "prob") %>%  select(prob, `2.5%`, `50%`, `97.5%`)

predict_data = predict %>% mutate(prob = str_replace_all(prob, "x_t_1\\[|\\]", "")) %>% 
  separate(prob, c("tp", "cluster"), ",") %>% mutate(time_points = as.numeric(tp)) %>% 
  left_join(new_time_index)

tot_growth = rstan::summary(run_fit, par = "y_t_out")$summary %>% as.tibble(rownames = "prob") %>%  
  select(prob, `2.5%`, `50%`, `97.5%`)  %>% mutate(time = new_time_points)

data_y_in = tibble(time = time_points, total = data_y) #, std = data_y_sd

pp2 = 
  
  tot_growth %>%
  ggplot(aes(x = time)) +
  geom_point(data = data_y_in, aes(x = time, y = log10(total), colour = "Experiment"), size = 3) + #colour = "red"
  #geom_errorbar(data = data_y_in, aes(x = time, ymax= total+1.9*std, ymin= total-1.9*std), width =0.2, colour = "red")+
  geom_line(aes(y = `50%`, linetype = "Modelling"), size = 0.65) + 
  geom_ribbon(aes(ymin= `2.5%`, ymax =`97.5%`, size = "Credible interval"), alpha = 0.1, fill ="blue") +
  geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  #facet_wrap(~cluster, nrow = 1) +
  #scale_y_log10()+
  theme_bw() +
  theme(text = element_text(size=15), legend.title=element_blank()) +
  ylab("Total Number of Cells") +
  xlab("time") 


time_index = data.frame( time = c(1:length(time_points)),  tp = time_points) %>% 
  mutate(tp = as.numeric(tp)) %>% as.tibble()

state_data = data_t %>% t() %>% as.tibble %>% mutate(cluster = 1:n()) %>% 
  gather(time, prob_val, -cluster) %>% mutate(time = str_replace_all(time, "V", ""), cluster = as.character(cluster)) %>% 
  mutate(time = as.integer(time)) %>% 
  left_join(time_index) %>%  mutate(cluster = as.integer(cluster))

#####################print the simulation data########################
predict = predict_data

state_data_t = state_data

cluster_name = tibble(cluster = 1:4, cluster_name = factor(c("CD34+", "CD34-", "CD235ab+", "HBA+"), levels = c("CD34+", "CD34-", "CD235ab+", "HBA+")))

# predict  %>% 
#   mutate(cluster = as.integer(cluster), Exp = factor(Exp, levels = c("BirA", "BioG1", "BioG1s"))) %>% #View()
#   select(-tp) %>% 
#   left_join(time_index) %>% 
#   left_join(cluster_name) %>% #View()
#   ggplot(aes(x = tp, y = `50%`)) +
#   geom_point(data = state_data_t %>% left_join(cluster_name) %>% mutate(Exp = factor(Exp, levels = c("BirA", "BioG1", "BioG1s"))), aes(x = tp, y = prob_val, colour = "Experiment")) + #colour = "red"
#   geom_line(aes(linetype = "Modelling"), size = 1.2) + 
#   geom_ribbon(aes(ymin= `2.5%`, ymax =`97.5%`, size = "Credible interval"), alpha = 0.1, fill ="blue") +
#   #geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
#   facet_grid(Exp~cluster_name) +
#   theme_light() +
#   theme(text = element_text(size=15), legend.title=element_blank()) +
#   ylab("Proportion of Cells") +
#   xlab("Day") 

pp3 = predict_data %>% 
  mutate(cluster = as.integer(cluster)) %>% 
  select(-tp) %>% 
  #left_join(time_index) %>% 
  left_join(cluster_name) %>% 
  ggplot(aes(x = time, y = `50%`)) +
  geom_point(data = state_data %>% left_join(cluster_name), aes(x = tp, y = prob_val, colour = "Experiment")) + #colour = "red"
  geom_line(aes(linetype = "Modelling"), size = 0.65) + 
  geom_ribbon(aes(ymin= `2.5%`, ymax =`97.5%`, size = "Credible interval"), alpha = 0.2, fill ="blue") +
  #geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  facet_wrap(~cluster_name, nrow = 1) +
  theme_classic() +
  theme(text = element_text(size=15), legend.title=element_blank()) +
  ylab("Proportion of Cells") +
  xlab("Day") 

gridExtra::grid.arrange(pp2, pp3,ncol = 1)

alpha_list1 = 
  alpha1$summary %>% as.tibble(rownames = "prob") %>% select(prob, `97.5%`) %>% 
  mutate(prob = str_replace_all(prob, "alpha1\\[|\\]", "")) %>% 
  separate(prob, c("To", "From"), ",") %>% 
  mutate(`50%` = if_else(To==From, 0, `97.5%`)) %>% 
  #mutate(`2.5%` = if_else(To==From, 0, `2.5%`)) %>% 
  mutate(To = as.numeric(To), From = as.numeric(From)) %>% 
  rename(lambda1 = `50%`) %>% View()

alpha_list1 = 
  alpha1$summary %>% as.tibble(rownames = "prob") %>% select(prob, `2.5%`, `50%`, `97.5%`) %>% 
  mutate(prob = str_replace_all(prob, "alpha1\\[|\\]", "")) %>% 
  separate(prob, c("To", "From"), ",") %>% #View()
mutate(`50%` = if_else(To==From, 0, `97.5%`)) %>% 
  #mutate(`2.5%` = if_else(To==From, 0, `2.5%`)) %>% 
  mutate(To = as.numeric(To), From = as.numeric(From)) %>% 
  rename(lambda1 = `50%`) %>% View()

lambda_list1 = 
  rstan::summary(run_fit, par = c("lambda1", "lambda_net"))$summary %>% as.tibble(rownames = "prob") %>% 
  select(prob, `2.5%`, `50%`, `97.5%`)

lambda_list2 = 
  rstan::summary(run_fit, par = c("lambda1", "lambda_net"))$summary %>% as.tibble(rownames = "prob") %>% 
  select(prob, `2.5%`, `50%`, `97.5%`)

lambda_list3 = 
  rstan::summary(run_fit, par = c("lambda1", "lambda_net"))$summary %>% as.tibble(rownames = "prob") %>% 
  select(prob, `2.5%`, `50%`, `97.5%`)

# alpha_list1 %>% 
#   ggplot(aes(From, To))+
#   geom_tile(aes(fill = lambda1), colour = "white")+
#   scale_fill_gradient(low = "white", high = "steelblue")

alpha1$summary %>% as.tibble(rownames = "prob") %>% View()

#########################################

##for the ICL
calculate_the_loglikeli2 <- function(c, r, test_data) {
  
  dir_ = paste("Results_by_date/Oct2019/full_protein_init_c", c, "_", r, ".rda", sep = "")
  run_fit = loadRData(dir_)
  
  vals = run_fit$par$log_like_i_gen %>% sum()
  
  return(vals)
  
}

calculate_the_loglikeli3 <- function(c, r, test_data) {
  
  dir_ = paste("Results_by_date/Oct2019/full_protein_init_c", c, "_", r, ".rda", sep = "")
  run_fit = loadRData(dir_)
  
  vals = run_fit$par$log_like_sum_gen %>% sum()
  
  return(vals)
  
}

calculate_the_loglikeli3(20, 1, sub_data)

no_cluster = c(10:25)
no_run = c(1:16)

log_prob = c()
best_log_prob = c()

for (k in no_cluster) {
  
  for(l in no_run) {
    
    p = calculate_the_loglikeli2(k, l, sub_data)
    
    best_log_prob = c(best_log_prob, p)
    
  }
}

for (k in no_cluster) {
  
  for(l in no_run) {
    
    p = calculate_the_loglikeli3(k, l, sub_data)
    
    log_prob = c(log_prob, p)
  }
}

# log_like_prob = tibble(no_cluster = c(1:6)) %>% mutate(log_prob = log_prob)

log_like_prob = log_prob %>% matrix(ncol = length(no_run), byrow = T) %>% as.tibble() %>% 
  mutate(no_cluster = no_cluster) %>% 
  gather(run, value, -no_cluster)

best_log_like_prob = best_log_prob %>% matrix(ncol = length(no_run), byrow = T) %>% as.tibble() %>% 
  mutate(no_cluster = no_cluster) %>% 
  gather(run, value, -no_cluster)

sample_size = log(length(sub_data))
day_vals = vec_d %>% unique()

log_like_prob %>% mutate(BIC_ = -2*value + 
                           sample_size* (no_cluster * nrow(prot_list) + no_cluster * nrow(prot_list) + no_cluster* length(day_vals))) %>% 
  
  mutate(ICL_ = - best_log_like_prob$value + 
           0.5* sample_size* (no_cluster * nrow(prot_list) + no_cluster * nrow(prot_list) + no_cluster* length(day_vals))) %>% 
  select(-value) %>% 
  gather(criteria, value, BIC_:ICL_) %>% 
  # group_by(no_cluster, criteria) %>% 
  # summarise(loglike = min(value)) %>% 
  ggplot(aes(x = no_cluster, y = value))+
  geom_point(aes(colour = criteria), size = 2.5) +
  #xlim(3.5, NA) +
  facet_wrap(~criteria, scales = "free") 

log_like_prob %>% mutate(BIC_ = -2*value + 
                           sample_size* (no_cluster * nrow(prot_list) + no_cluster * nrow(prot_list) + no_cluster* length(day_vals))) %>% 
  
  mutate(ICL_ = - best_log_like_prob$value + 
           0.5* sample_size* (no_cluster * nrow(prot_list) + no_cluster * nrow(prot_list) + no_cluster* length(day_vals))) %>% 
  select(-value) %>% 
  gather(criteria, value, BIC_:ICL_) %>% 
  group_by(no_cluster, criteria) %>%
  summarise(loglike = min(value)) %>%
  ggplot(aes(x = no_cluster, y = loglike))+
  geom_point(aes(colour = criteria), size = 2.5) +
  #xlim(3.5, NA) +
  facet_wrap(~criteria, scales = "free") 

no_cluster = 25
no_run = c(1:16)

log_prob = c()
best_log_prob = c()

for (k in no_cluster) {
  
  for(l in no_run) {
    
    p = calculate_the_loglikeli2(k, l, sub_data)
    
    best_log_prob = c(best_log_prob, p)
    
  }
}

for (k in no_cluster) {
  
  for(l in no_run) {
    
    p = calculate_the_loglikeli3(k, l, sub_data)
    
    log_prob = c(log_prob, p)
  }
}

best_log_prob %>% as.tibble() %>% mutate(no_run) %>% 
  arrange(value)
