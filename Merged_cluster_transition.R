library(rstan)
library(tidyverse)

############find the knn list#########
cluster_prob = all_cells_prob 

corre_data = matrix(0, ncol = 5, nrow =ncol(cluster_prob)*ncol(cluster_prob))

for(k in 1:ncol(cluster_prob)){ #!k %in% empty_cluster
  
  for(j in 1:ncol(cluster_prob)) {
    
    corr =0
    count_p1 = 0
    count_p2 = 0
    
    for(i in 1:nrow(cluster_prob)){
      
      p1 = cluster_prob[i, k]
      p2 = cluster_prob[i, j]
      pp = p1 + p2
      
      if(pp>=0.95) {
        
        if(p1 > 0.95) count_p1 = count_p1 + 1
        else if(p2 > 0.95) count_p2 = count_p2 + 1
        else corr = corr + 1
        
      }
    }
    
    index = ncol(cluster_prob) * (k-1) + j
    
    print(index)
    
    sum = count_p1 + count_p2 + corr
    
    corre_data[index, 1] =  k
    corre_data[index, 2] =  j
    corre_data[index, 3] =  count_p1/sum
    corre_data[index, 4] =  count_p2/sum
    corre_data[index, 5] =  corr/sum
  }
}

knn_list = corre_data %>% as.tibble() %>% rename(cluster = V1, cluster_n = V2, p1 = V3, p2 = V4, p = V5) %>%
  select(-p1, -p2) %>% group_by(cluster) %>% arrange(cluster, desc(p)) %>% ungroup() %>% filter(p>=1e-2) %>% 
  select(cluster, cluster_n) %>%
  filter(cluster != cluster_n) %>%
  rename(state.x = cluster, state.y = cluster_n)

############CALCULATE CORRELATIONS between clusters##############

cluster = c(1:25)

corr_mean = c()
corr_prop = c()

x_mean = signal %>% select(cluster, prot_name, vals) %>% 
  filter(!prot_name %in% c("IdU", "mean_Ir", "p_Rb"))

x_prop = data_x[-650918,] %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
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
  xlim(0.9, 1.0)+
  ylim(0.9, 1.0)

corr_mean %>% mutate(corr_mu = corr, corr_prop = corr_prop$corr) %>%
  mutate(cluster = paste0(cluster1, " ", cluster2)) %>% 
  filter(cluster1 != cluster2) %>% 
  filter(cluster1 <= cluster2) %>% 
  filter(corr_mu >=0.9, corr_prop >= 0.9) %>% 
  select(cluster1, cluster2) %>% 
  rename(cluster = cluster2, merge_cluster = cluster1)

##################################C=25
cluster_merge = 
  tibble(cluster = c(1:25),
         merge_cluster = c(1, NA, 3, 4, 5, 6, 7, 8, 5, 6, 11, 5, 1,
                           NA, 15, 16, 17, 18, 4, 3, 17, 4, 11, 24, 18))

merge_cluster = c(1, 3, 4, 5, 6, 7, 8, 5, 6, 11, 5, 1,
                  15, 16, 17, 18, 4, 3, 17, 4, 11, 24, 18) %>% unique()
cluster_new = 
  tibble(merge_cluster = merge_cluster,
         state = c(1:length(merge_cluster)))

##################################C=30
cluster_merge = 
  tibble(cluster = c(1:30),
         merge_cluster = c(1, 1, 3, 4, 5, 6, NA, 8, 9, 10, 1, NA, NA, 14, NA,
                           9, 17, 2, 19, 9, 4, 19, 23, 24, NA, 3, 27, 28, 14, 9))

merge_cluster = c(1, 1, 3, 4, 5, 6, 8, 9, 10, 1, 14, 9, 17, 2, 19, 9, 4, 19, 23, 24, 3, 27, 28, 14, 9) %>% unique()

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

p2 = data_x[-650918,] %>% mutate(cluster = as.integer(all_cells)) %>%
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

###################modify the date with new clusters###############################
data_t =  data_x[-650918,] %>% mutate(cluster = all_cells) %>% 
  left_join(cluster_merge %>% mutate(cluster = as.numeric(cluster)), by = c("cluster" = "cluster")) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>%
  group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
  filter(Day != 0) %>% #Day != 11
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(state, prob, fill = 0) %>%
  select(-Day) %>% 
  as.matrix()

knn_list = knn_list %>% 
  left_join(cluster_merge, by = c("state.x" = "cluster")) %>% 
  left_join(cluster_merge, by = c("state.y" = "cluster")) %>% 
  filter(merge_cluster.x != "NA", merge_cluster.y != "NA") %>% 
  select(-state.x, -state.y) %>% 
  left_join(cluster_new, by = c("merge_cluster.x" = "merge_cluster")) %>%
  left_join(cluster_new, by = c("merge_cluster.y" = "merge_cluster")) %>% 
  select(state.x, state.y) %>% 
  filter(state.x != state.y) %>% 
  unique()

n_state = nrow(cluster_new)
n_time = nrow(data_t)
n_rate = nrow(knn_list)
n_trans = 3

time_points = data_x[-650918,] %>% select(Day) %>% filter(Day !=0) %>% unique() %>% pull(Day) %>% as.integer()

growth_total = c(14.7, 70.42, 295.76, 921.1, 5205.06, 5619.63, 9553.36, 14610.96, 11044.95, 8242.5, 6594.) * 1e6 #7 at Day0

data_y = growth_total %>% ceiling()

typeof(data_y)

#########S phase ratio calculation################

S_Phase_ratio = data_x[-650918,] %>% mutate(cluster = all_cells) %>% 
  left_join(cluster_merge %>% mutate(cluster = as.numeric(cluster)), by = c("cluster" = "cluster")) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>%
  select(Day, IdU, state) %>%
  mutate(S_Phase = ifelse(IdU <= sinh(3.466711), "Other", "S_Phase")) %>% 
  group_by(state) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  group_by(state, S_Phase) %>% 
  mutate(count = n()/total) %>% ungroup() %>% 
  select(state, S_Phase, count) %>% 
  unique() %>% 
  filter(S_Phase == "S_Phase") %>% 
  arrange(state) %>% 
  select(count) %>%
  pull(count)


p3 = data_x[-650918,] %>% mutate(cluster = all_cells) %>% 
  left_join(cluster_merge %>% mutate(cluster = as.numeric(cluster)), by = c("cluster" = "cluster")) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>%
  select(Day, IdU, state) %>%
  mutate(S_Phase = ifelse(IdU <= sinh(3.5), "Other", "S_Phase")) %>% 
  group_by(Day, state) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  group_by(state, Day, S_Phase) %>% 
  mutate(count = n()/total) %>% ungroup() %>% 
  select(Day, state, S_Phase, count) %>% 
  unique() %>% 
  filter(S_Phase == "S_Phase") %>%
  spread(Day, count, fill = 0) %>% 
  gather(Day, count, -state, -S_Phase) %>% 
  mutate(Day = as.numeric(Day)) %>% 
  ggplot()+
  geom_point(aes(x= Day, y = count), colour = "blue", size = 2.5) +
  facet_wrap(~state, nrow = 2) +
  theme_bw()

gridExtra::grid.arrange(p2, p3, ncol =1)

S_Phase_ratio = data_x[-650918,] %>% mutate(cluster = all_cells) %>% 
  left_join(cluster_merge %>% mutate(cluster = as.numeric(cluster)), by = c("cluster" = "cluster")) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>%
  select(Day, IdU, state) %>%
  mutate(S_Phase = ifelse(IdU <= sinh(3.5), "Other", "S_Phase")) %>%
  mutate(
    Time_inter = case_when(
      Day < 11   ~  1,
      Day >= 11 & Day < 15 ~ 2,
      Day >= 15 ~  3
    )
  ) %>% 
  group_by(Time_inter, state) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  group_by(Time_inter, state, S_Phase) %>% 
  mutate(count = n()/total) %>% ungroup() %>% 
  select(state, Time_inter, S_Phase, count) %>%
  unique() %>% 
  filter(S_Phase == "S_Phase") %>% 
  spread(Time_inter, count, fill = 0) %>% 
  gather(Time_inter, count, -state, -S_Phase) %>% 
  arrange(Time_inter, state) %>%
  pull(count)
############dirichlet fit############  
cmp_model = stan_model(file = "Results_by_date/Aug2019/State_transitions_w_NetGrowth_test_v6.stan") #Results_by_date/Aug2019/State_transitions_w_NetGrowth.stan
cmp_model = stan_model(file = "Software/cluster_model/State_transitions_w_NetGrowth.stan") #Results_by_date/Aug2019/State_transitions_w_NetGrowth.stan


gen_init_vals = function(n_rate, n_state, n_time){
  list(#lambda1 = runif(n = n_rate, 1, 2),
    #lambda_death = runif(n = 3*n_state, 0, 0),
    Sigma_vec1 = runif(n = n_state, 1000, 3000),
    Sigma_vec2 = runif(n = n_state, 1000, 3000),
    Sigma_vec3 = runif(n = n_state, 1000, 3000),
    lambda_trans1 = runif(n = 1, 1.5, 2),
    lambda_trans2 = runif(n = 1, 1.5, 2),
    Sigma_tot = runif(n = 1, 0, 1))
} # true_y = runif(n = n_time, 1e5, 5e5)

fit_Yule_func = function(cmp_model, data_t, data_y, n_state, n_time, time_points, n_rate, knn_list, n_trans, num_chains = 4, iter = 500){
  
  n_state = ncol(data_t)
  n_time = nrow(data_t)
  n_rate = nrow(knn_list)
  
  data_for_stan = list(n_state = n_state, n_time = n_time, x_t_in = data_t, y_t_in = data_y, #y_t_sd_in = data_y_sd, 
                       time_points = time_points, n_rate_in=n_rate, knn_weight = knn_list,
                       n_trans = n_trans)
                       #s_phase_ratio = S_Phase_ratio) 
  
  init_vals = rerun(num_chains, gen_init_vals(n_rate, n_state, n_time))
  
  run_fit = sampling(cmp_model, iter = iter, data = data_for_stan, chains = num_chains, thin = 1, cores = 14, init = init_vals) # warmup = 500, init = init_vals
  
  run_fit
}

run_fit = fit_Yule_func(cmp_model, data_t, data_y, n_state, n_time, time_points, n_rate, knn_list, n_trans)

save(run_fit, file = "Results_by_date/Dec2019/lambda_trans_with_gen_quant_13states_new.rda")
load(file = "Results_by_date/Dec2019/3_lambda_trans.rda")

traceplot(run_fit, inc_warmup = F, pars = "lambda_death")

p1 = plot(run_fit, pars = "lambda1")
p2 = plot(run_fit, pars = "lambda2")
p3 = plot(run_fit, pars = "lambda3")

tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
tbl_summary %>% View()
alpha1 = rstan::summary(run_fit, par = "alpha1")
alpha2 = rstan::summary(run_fit, par = "alpha2")
alpha3 = rstan::summary(run_fit, par = "alpha3")


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

pp1 = tot_growth %>%
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

predict = predict_data

state_data_t = state_data

cluster_name = tibble(cluster = 1:4, cluster_name = factor(c("CD34+", "CD34-", "CD235ab+", "HBA+"), levels = c("CD34+", "CD34-", "CD235ab+", "HBA+")))

pp2 = predict_data %>% 
  mutate(cluster = as.integer(cluster)) %>% 
  select(-tp) %>% 
  #left_join(time_index) %>%
  ggplot(aes(x = time, y = `50%`)) +
  geom_point(data = state_data, aes(x = tp, y = prob_val, colour = "Experiment")) + #colour = "red"
  geom_line(aes(linetype = "Modelling"), size = 0.65) + 
  geom_ribbon(aes(ymin= `2.5%`, ymax =`97.5%`, size = "Credible interval"), alpha = 0.2, fill ="blue") +
  #geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  facet_wrap(~cluster, nrow = 2) +
  theme_classic() +
  theme(text = element_text(size=15), legend.title=element_blank()) +
  ylab("Proportion of Cells") +
  xlab("Day") 

gridExtra::grid.arrange(pp1, pp2,ncol = 1)

# data1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>%
#   left_join(cluster_merge) %>% 
#   filter(merge_cluster != "NA") %>% 
#   left_join(cluster_new) 
# 
# data2 = data_x %>% mutate(cluster = as.integer(all_cells)) %>%
#   left_join(cluster_merge) %>% 
#   filter(merge_cluster != "NA") %>% 
#   left_join(cluster_new) %>% 
#   group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
#   mutate(prob = num/sum(num)) %>% ungroup() %>% 
#   select(-num) %>% 
#   spread(state, prob, fill = 0) %>% 
#   gather(state, prob, -Day) %>% 
#   mutate(state = as.numeric(state))
# 
# save(data1, file = "Results_by_date/Dec2019/prot_data.rda")
# save(data2, file = "Results_by_date/Dec2019/prop_data.rda")
# 
# load("Results_by_date/Dec2019/prot_data.rda")
# load("Results_by_date/Dec2019/prop_data.rda")

####S_ratio model

cmp_model = stan_model(file = "Results_by_date/Aug2019/State_transitions_w_NetGrowth_test_v6.stan")

fit_Yule_func = function(cmp_model, data_t, data_y, n_state, n_time, time_points, n_rate, knn_list, S_Phase_ratio, n_trans, num_chains = 4, iter = 100){
  
  n_state = ncol(data_t)
  n_time = nrow(data_t)
  n_rate = nrow(knn_list)
  
  data_for_stan = list(n_state = n_state, n_time = n_time, x_t_in = data_t, y_t_in = data_y, #y_t_sd_in = data_y_sd, 
                       time_points = time_points, n_rate_in=n_rate, knn_weight = knn_list,
                       n_trans = n_trans, s_phase_ratio = S_Phase_ratio)
  #s_phase_ratio = S_Phase_ratio) 
  
  init_vals = rerun(num_chains, gen_init_vals(n_rate, n_state, n_time))
  
  run_fit = sampling(cmp_model, iter = iter, data = data_for_stan, chains = num_chains, thin = 1, cores = 14, init = init_vals) # warmup = 500, init = init_vals
  
  run_fit
}

run_fit = fit_Yule_func(cmp_model, data_t, data_y, n_state, n_time, time_points, n_rate, knn_list, S_Phase_ratio, n_trans)
 
save(run_fit, file = "Results_by_date/Results2020/lambda_trans_with_gen_quant_0213.rda")
load(file = "Results_by_date/Dec2019/lambda_trans_with_gen_quant_13states_13*3_no_S_50knn.rda")

traceplot(run_fit, inc_warmup = F, pars = "lambda_death")
plot(run_fit, pars = "lambda_death")

p1 = plot(run_fit, pars = "lambda1")
p2 = plot(run_fit, pars = "lambda2")
p3 = plot(run_fit, pars = "lambda3")

tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
tbl_summary %>% View()

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

pp1 = tot_growth %>%
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

predict = predict_data

state_data_t = state_data

pp2 = predict_data %>% 
  mutate(cluster = as.integer(cluster)) %>% 
  select(-tp) %>% 
  #left_join(time_index) %>%
  ggplot(aes(x = time, y = `50%`)) +
  geom_point(data = state_data, aes(x = tp, y = prob_val, colour = "Experiment")) + #colour = "red"
  geom_line(aes(linetype = "Modelling"), size = 0.65) + 
  geom_ribbon(aes(ymin= `2.5%`, ymax =`97.5%`, size = "Credible interval"), alpha = 0.2, fill ="blue") +
  #geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  facet_wrap(~cluster, nrow = 2) +
  theme_classic() +
  theme(text = element_text(size=15), legend.title=element_blank()) +
  ylab("Proportion of Cells") +
  xlab("Day") 

gridExtra::grid.arrange(pp1, pp2,ncol = 1)

vals = map_est4$par 

theta = tibble(vals = vals, param = names(vals)) %>% 
  mutate(param = str_remove_all(param, "\\]")) %>% separate(param, c("var", "c_p"), sep = "\\[") %>% 
  separate(c_p, c("tp_cluster", "protein"), sep = "\\,") 

lambda1 = rstan::summary(run_fit, par = "lambda1")
lambda2 = rstan::summary(run_fit, par = "lambda2")
lambda3 = rstan::summary(run_fit, par = "lambda3")

################knn network plot #######################
library(RANN)
library(tidygraph)
library(ggraph)

lambda_diff1 = data_x[-650918,] %>% mutate(cluster = all_cells) %>% 
  left_join(cluster_merge %>% mutate(cluster = as.numeric(cluster)), by = c("cluster" = "cluster")) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>%
  group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(Day, prob, fill = 0) %>% 
  select(`2`:`11`, state) %>% 
  gather(Day, prob, -state) %>% 
  group_by(state) %>% 
  mutate(diff = max(prob) - min(prob)) %>% 
  ungroup() %>% 
  select(state, diff) %>% 
  unique()

lambda_diff2 = data_x[-650918,] %>% mutate(cluster = all_cells) %>% 
  left_join(cluster_merge %>% mutate(cluster = as.numeric(cluster)), by = c("cluster" = "cluster")) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>%
  group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(Day, prob, fill = 0) %>% 
  select(`11`:`14`, state) %>% 
  gather(Day, prob, -state) %>% 
  group_by(state) %>% 
  mutate(diff = max(prob) - min(prob)) %>% 
  ungroup() %>% 
  select(state, diff) %>% 
  unique()
  
lambda_diff3 = data_x[-650918,] %>% mutate(cluster = all_cells) %>% 
  left_join(cluster_merge %>% mutate(cluster = as.numeric(cluster)), by = c("cluster" = "cluster")) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>%
  group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(Day, prob, fill = 0) %>% 
  select(`16`:`20`, state) %>% 
  gather(Day, prob, -state) %>% 
  group_by(state) %>% 
  mutate(diff = max(prob) - min(prob)) %>% 
  ungroup() %>% 
  select(state, diff) %>% 
  unique()

lambda_diff = lambda_diff1 %>% rename(diff1 = diff) %>% mutate(diff2 = lambda_diff2$diff, diff3 = lambda_diff3$diff)

####death or birth rate
growth_net = rstan::summary(run_fit, par = "lambda_death") 

#growth_net = #
growth_net$summary %>% as.tibble(rownames = "rate") %>% 
  pull(`50%`) %>% 
  matrix(ncol = 3) %>% 
  as.tibble() %>% 
  mutate(state = 1:c) %>% 
  rename(t1 = V1, t2 = V2, t3 = V3) %>% 
  
#find the transition rate matrix
alpha_list1 = lambda1$summary %>% as.tibble(rownames = "prob") %>% 
  select(`50%`) %>% 
  cbind(knn_list) %>% 
  rename(trans_rate = `50%`, To = state.x, From = state.y) %>% 
  mutate(To = as.numeric(To), From = as.numeric(From)) %>% 
  as.tibble()

alpha_list2 = lambda2$summary %>% as.tibble(rownames = "prob") %>% 
  select(`50%`) %>% 
  cbind(knn_list) %>% 
  rename(trans_rate = `50%`, To = state.x, From = state.y) %>% 
  mutate(To = as.numeric(To), From = as.numeric(From)) %>% 
  as.tibble()

alpha_list3 = lambda3$summary %>% as.tibble(rownames = "prob") %>% 
  select(`50%`) %>% 
  cbind(knn_list) %>% 
  rename(trans_rate = `50%`, To = state.x, From = state.y) %>% 
  mutate(To = as.numeric(To), From = as.numeric(From)) %>% 
  as.tibble()

connect_list = knn_list

######### draw the full network ##############    
network1 = alpha_list1 %>% 
  left_join(lambda_diff, by = c("From" = "state")) %>% 
  rename(diff_From1 = diff1) %>% 
  left_join(lambda_diff, by = c("To" = "state")) %>% 
  rename(diff_To1 = diff1) %>% 
  mutate(new_lambda1 = trans_rate * diff_From1 / (diff_To1+1e-10)) %>% 
  mutate(new_lambda1 = if_else(diff_To1 <= 1e-2 | diff_From1 <= 1e-2, 0, new_lambda1)) %>% 
  filter(new_lambda1 > 0.01) %>% select(From, To, new_lambda1)

network2 = alpha_list2 %>% 
  left_join(lambda_diff, by = c("From" = "state")) %>% 
  rename(diff_From1 = diff2) %>% 
  left_join(lambda_diff, by = c("To" = "state")) %>% 
  rename(diff_To1 = diff2) %>% 
  mutate(new_lambda1 = trans_rate * diff_From1 / (diff_To1+1e-10)) %>% 
  mutate(new_lambda1 = if_else(diff_To1 <= 1e-2 | diff_From1 <= 1e-2, 0, new_lambda1)) %>% 
  filter(new_lambda1 > 0.01) %>% select(From, To, new_lambda1)

network3 = alpha_list3 %>% 
  left_join(lambda_diff, by = c("From" = "state")) %>% 
  rename(diff_From1 = diff3) %>% 
  left_join(lambda_diff, by = c("To" = "state")) %>% 
  rename(diff_To1 = diff3) %>% 
  mutate(new_lambda1 = trans_rate * diff_From1 / (diff_To1+1e-10)) %>% 
  mutate(new_lambda1 = if_else(diff_To1 <= 1e-2 | diff_From1 <= 1e-2, 0, new_lambda1)) %>% 
  #filter(new_lambda1 > 0.01) %>% 
  select(From, To, new_lambda1)

c = 13

###############from 2 to 11###########

cluster_data = data_x[-650918,] %>% mutate(cluster = as.integer(all_cells)) %>%
  left_join(cluster_merge) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>% 
  group_by(Day, state) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(state, prob, fill = 0) %>% 
  gather(state, prob, -Day) %>% 
  mutate(state = as.numeric(state))
  
  
useful_cluster = cluster_data  %>%
  filter(Day >= 15) %>%  #Day >0, Day <= 11
  filter(prob> 0) %>% 
  select(state) %>%
  filter(state != 6) %>% unique() %>% pull()
  

alpha_list3 %>%
  left_join(network3, by = c("From" = "From", "To" = "To")) %>% 
  replace_na(replace = list(new_lambda1 = 0)) %>% 
  # mutate(new_lambda = (new_lambda1 + new_lambda2)/2) %>% 
  # #mutate(new_lambda = if_else(new_lambda1 == 0, 0, new_lambda1 + new_lambda2)) %>% 
  # #mutate(new_lambda = if_else(new_lambda2 == 0, 0, new_lambda1 + new_lambda2)) %>% 
  filter(new_lambda1 != 0) %>%
  #mutate(new_lambda1 = 100*new_lambda1) %>% 
  select(From, To, trans_rate, new_lambda1) %>%
  filter(From %in% useful_cluster, To %in% useful_cluster) %>% 
  as_tbl_graph() %>% 
  activate(nodes) %>% 
  mutate(name = as.integer(name)) %>% 
  arrange(name) %>% 
  ggraph(layout = "auto") +#, node.positions = node_pos) +
  geom_edge_link(aes(alpha = trans_rate), arrow = arrow(length = unit(7.5, 'mm')), colour = "blue") +
  # geom_edge_link(aes(alpha = new_lambda2), arrow = arrow(length = unit(5.5, 'mm')), colour = "blue") +
  geom_node_point(aes(col = factor(name, levels = 1:c)), size = 10 ) + geom_node_text(aes(label = name), size = 7.5) +
  #theme(panel.background = element_rect(color = "grey"))
  theme_minimal()+
  theme(legend.position = "none")

network3 %>%
  left_join(alpha_list1, by = c("From" = "From", "To" = "To")) %>% 
  left_join(alpha_list2, by = c("From" = "From", "To" = "To")) %>% 
  left_join(alpha_list3, by = c("From" = "From", "To" = "To")) %>% 
  left_join(network1, by = c("From" = "From", "To" = "To")) %>% 
  left_join(network2, by = c("From" = "From", "To" = "To")) %>% 
  replace_na(replace = list(new_lambda1.x = 0)) %>% 
  replace_na(replace = list(new_lambda1.y = 0)) %>% 
  replace_na(replace = list(new_lambda1 = 0)) %>% 
  mutate(total_lambda = new_lambda1.x + new_lambda1.y + new_lambda1) %>% 
  filter(total_lambda != 0) %>%
  select(From, To, trans_rate, trans_rate.x, trans_rate.y) %>%
  filter(From %in% useful_cluster, To %in% useful_cluster) %>% 
  as_tbl_graph() %>% 
  activate(nodes) %>% 
  mutate(name = as.integer(name)) %>% 
  arrange(name) %>% 
  ggraph(layout = "auto") +#, node.positions = node_pos) +
  geom_edge_link(aes(alpha = trans_rate), arrow = arrow(length = unit(7.5, 'mm')), colour = "blue") +
  geom_edge_link(aes(alpha = trans_rate.x), arrow = arrow(length = unit(7.5, 'mm')), colour = "red") +
  geom_edge_link(aes(alpha = trans_rate.y), arrow = arrow(length = unit(7.5, 'mm')), colour = "green") +
  geom_node_point(aes(col = factor(name, levels = 1:c)), size = 10 ) + geom_node_text(aes(label = name), size = 7.5) +
  theme_minimal()+
  theme(legend.position = "none")

############Plot the heatmap for clusters#########

signal %>%
  left_join(cluster_merge) %>% 
  filter(merge_cluster != "NA") %>% 
  left_join(cluster_new) %>% 
  select(state, prot_name, vals) %>% 
  group_by(state, prot_name) %>% 
  mutate(value = mean(vals)) %>% 
  select(state, prot_name, value) %>% 
  unique() %>% 
  ungroup() %>%
  spread(state, value) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "prot_name") %>% 
  asinh() %>% pheatmap::pheatmap()
 
node_pos = matrix(c(1, 0, 10, 
                    2, 0, 5,
                    3, 0, 0,
                    4, 5, 0,
                    5, 5, 5,
                    6, 5, 10,
                    7, 10, -2.5,
                    8, 10, 2.5,
                    9, 10, 10,
                    10, 10, 5,
                    11, 15, 5,
                    12, 15, 10,
                    13, 17.5, 5,
                    14, 17.5, 10,
                    16, 20, 5,
                    17, 22.5, 5,
                    18, 25, 5), ncol = 3, byrow = T) %>% 
  as.tibble() %>% rename(name = V1, x = V2, y = V3) %>%
  #mutate(name = as.character(name)) %>% 
  arrange(name) %>% 
  select(x, y)
###############################################
knn_list %>% as_tbl_graph() %>% 
  #activate(nodes) %>% 
  #left_join(node_pos) %>% 
  ggraph(layout = "auto") +
  geom_edge_link(arrow = arrow(length = unit(5, 'mm'))) +
  geom_node_point(aes(col = factor(name, levels = 1:c)), size = 8) + geom_node_text(aes(label = name))


## knn network all genes
tdy_knn = tidy_knn(full_connect_list)

#plot knn
ggraph(tdy_knn, layout = "fr") +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(col = name, size = 0.1))

connect_list %>% 
  left_join(alpha_list2) %>% 
  left_join(lambda_diff, by = c("From" = "cluster")) %>% 
  rename(diff_From = diff2) %>% 
  left_join(lambda_diff, by = c("To" = "cluster")) %>% 
  rename(diff_To = diff2) %>% 
  mutate(new_lambda = lambda2 * diff_From / (diff_To+1e-10)) %>% 
  mutate(new_lambda = if_else(diff_To <= 1e-2 | diff_From <= 1e-2, 0, new_lambda)) %>% #View()
  filter(new_lambda > 0.01) %>% 
  filter(!To %in% c(13, 20, 15), From != 13, From != 15) %>% 
  as_tbl_graph() %>% 
  #activate(nodes) %>% 
  #left_join(node_pos) %>% 
  ggraph(layout = "auto") +
  geom_edge_link(aes(alpha = new_lambda), arrow = arrow(length = unit(5, 'mm'))) +
  geom_node_point(aes(col = factor(name, levels = 1:c)), size = 8) + geom_node_text(aes(label = name))

