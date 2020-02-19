//Prediction for Oct K1 data
//22/01/2020
//State transition Dirichlet distribution
//fit both cell sum and proportion data
//Adding ratio of S_Phase cells in

functions {

  matrix lambda_to_alpha(int[,] knn_weight, int n_state, vector lambda, vector lambda_net) {
    
    int len_knn;
    matrix[n_state, n_state] alpha; 

    len_knn = dims(knn_weight)[1];

    alpha = rep_matrix(0., n_state, n_state);
  
    for(j in 1:len_knn){
         int pos1 = knn_weight[j,1];
         int pos2 = knn_weight[j,2];
         
         alpha[pos1, pos2] = lambda[j];
     }
     
     for(i in 1:n_state) alpha[i,i] = 0; //lambda_net[i] - sum(alpha[1:n_state, i]);
     
    return(alpha); 
  }
  
  matrix lambda_to_alpha_3(int[,] knn_weight, int n_state, matrix alpha1, matrix alpha2, matrix alpha3, vector lambda_death, real lambda_trans1, real lambda_trans2) {

  int len_knn = 3;

  matrix[(len_knn*n_state), (len_knn*n_state)] alpha;

  alpha = rep_matrix(0., (len_knn*n_state), (len_knn*n_state));

  for(j in 1:n_state){
    for(i in 1:n_state){
      alpha[j, i] = alpha1[j, i];
    }
  }

  for(j in (n_state+1):(2*n_state)){
    for(i in (n_state+1):(2*n_state)){
      alpha[j, i] = alpha2[j-n_state, i-n_state];
    }
  }

  for(j in (2*n_state+1):(3*n_state)){
    for(i in (2*n_state+1):(3*n_state)){
      alpha[j, i] = alpha3[j-(2*n_state), i-(2*n_state)];
    }
  }


  for(j in 1:n_state){
    alpha[j+n_state, j] = lambda_trans1;
  }

  for(j in (1+n_state):(2*n_state)){
    alpha[j+n_state, j] = lambda_trans2;
  }

  //cell cycle - cell death - transitions
  for(i in 1:n_state) alpha[i,i] = lambda_death[i] - sum(alpha[1:(len_knn*n_state), i]);
  for(i in (1+n_state):(2*n_state)) alpha[i,i] = lambda_death[i] - sum(alpha[1:(len_knn*n_state), i]); //-n_state
  for(i in (1+2*n_state):(3*n_state)) alpha[i,i] = lambda_death[i] - sum(alpha[1:(len_knn*n_state), i]); //-(2*n_state)

  return(alpha);
}

  matrix lambda_to_alpha_4(int[,] knn_weight, int n_state, matrix alpha1, matrix alpha2, matrix alpha3, vector s_rate, vector lambda_death, real lambda_trans1, real lambda_trans2, vector s_phase_ratio) {

    int len_knn = 3;//dims(knn_weight)[1];
    // int indx_start;
    // int indx_end;
    //len_knn ;
    matrix[(len_knn*n_state), (len_knn*n_state)] alpha;

    alpha = rep_matrix(0., (len_knn*n_state), (len_knn*n_state));

    for(j in 1:n_state){
      for(i in 1:n_state){
        alpha[j, i] = alpha1[j, i];
    }
  }

  for(j in (n_state+1):(2*n_state)){
    for(i in (n_state+1):(2*n_state)){
      alpha[j, i] = alpha2[j-n_state, i-n_state];
    }
  }

   for(j in (2*n_state+1):(3*n_state)){
      for(i in (2*n_state+1):(3*n_state)){
        alpha[j, i] = alpha3[j-(2*n_state), i-(2*n_state)];
      }
   }

  for(j in 1:n_state){
    alpha[j+n_state, j] = lambda_trans1;
  }

  for(j in (1+n_state):(2*n_state)){
    alpha[j+n_state, j] = lambda_trans2;
  }

   //cell cycle - cell death - transitions
  for(i in 1:n_state) alpha[i,i] =  s_rate[i]*s_phase_ratio[i] - lambda_death[i]  - sum(alpha[1:(len_knn*n_state), i]);
  for(i in (1+n_state):(2*n_state)) alpha[i,i] =  s_rate[i]*s_phase_ratio[i] - lambda_death[i]  - sum(alpha[1:(len_knn*n_state), i]);
  for(i in (1+2*n_state):(3*n_state)) alpha[i,i] =  s_rate[i]*s_phase_ratio[i] - lambda_death[i] - sum(alpha[1:(len_knn*n_state), i]); //-(2*n_state)


  // indx_start = len_knn+1;
  // indx_end   = num_elements(lambda);

    return(alpha);
  }

  matrix lambda_to_alpha_5(int[,] knn_weight, int n_state, matrix alpha1, matrix alpha2, matrix alpha3, real s_rate, vector lambda_death, real lambda_trans1, real lambda_trans2, vector s_phase_ratio) {

    int len_knn = 3;//dims(knn_weight)[1];
    // int indx_start;
    // int indx_end;
    //len_knn ;
    matrix[(len_knn*n_state), (len_knn*n_state)] alpha;

    alpha = rep_matrix(0., (len_knn*n_state), (len_knn*n_state));

    for(j in 1:n_state){
      for(i in 1:n_state){
        alpha[j, i] = alpha1[j, i];
    }
  }

  for(j in (n_state+1):(2*n_state)){
    for(i in (n_state+1):(2*n_state)){
      alpha[j, i] = alpha2[j-n_state, i-n_state];
    }
  }

   for(j in (2*n_state+1):(3*n_state)){
      for(i in (2*n_state+1):(3*n_state)){
        alpha[j, i] = alpha3[j-(2*n_state), i-(2*n_state)];
      }
   }

  for(j in 1:n_state){
    alpha[j+n_state, j] = lambda_trans1;
  }

  for(j in (1+n_state):(2*n_state)){
    alpha[j+n_state, j] = lambda_trans2;
  }

   //cell cycle - cell death - transitions
  for(i in 1:n_state) alpha[i,i] =  s_rate*s_phase_ratio[i] - lambda_death[i]  - sum(alpha[1:(len_knn*n_state), i]);
  for(i in (1+n_state):(2*n_state)) alpha[i,i] =  s_rate*s_phase_ratio[i] - lambda_death[i]  - sum(alpha[1:(len_knn*n_state), i]);
  for(i in (1+2*n_state):(3*n_state)) alpha[i,i] =  s_rate*s_phase_ratio[i] - lambda_death[i] - sum(alpha[1:(len_knn*n_state), i]); //-(2*n_state)


  // indx_start = len_knn+1;
  // indx_end   = num_elements(lambda);

    return(alpha);
  }

}

data {

  int<lower=0> n_state; //number of states
  int<lower=0> n_time; //number of time points
  vector[n_state] x_t_in[n_time]; //data at t
  real y_t_in[n_time]; //total number of viable cells per time point
  //vector[n_time] y_t_sd_in; //total number of viable cells per time point
  vector[n_time] time_points; //time measured
  int<lower=0> n_rate_in;
  int<lower=0> n_trans;
  int knn_weight[n_rate_in, 2]; //knn weight for clusters
  vector[3*n_state] s_phase_ratio;
}

// make knn symmetric and make list yes
// make alpha transformed param yes
// lambda of length list knn yes 
// make jump point in one yes
// Compare theta to predicted number of cells in each cluster 
transformed data {
  vector[n_state] x_t[n_time]; //data at t
  vector[n_time] y_t;
  
  int n_rate;
  
  for(t in 1:n_time) {
    
   x_t[t] = x_t_in[t] + 1e-6;
   x_t[t] = x_t[t]/sum(x_t[t]);
   
   y_t[t] = log10(y_t_in[t]+1e-6);
  }
  
  n_rate = n_rate_in;

}

parameters {
  vector<lower = 0>[n_rate] lambda1;
  vector<lower = 0>[n_rate] lambda2;
  vector<lower = 0>[n_rate] lambda3;
  
  vector<lower = 0, upper = 3>[3*n_state] lambda_death;
  // vector<lower = -2, upper = 4>[3*n_state] lambda_death;
  //vector<lower = 0>[n_state] lambda_death;
  // vector<lower = 2.4, upper= 4.8>[n_state] s_rate;
  //vector<lower = 0, upper= 4.8>[3*n_state] s_rate;
  real<lower = 1, upper= 4.8> s_rate;
   
  real<lower = 0, upper = 3> lambda_trans1;
  real<lower = 0, upper = 2> lambda_trans2;
  
  // vector<lower = 30>[n_state] Sigma_vec;
  vector<lower = 30>[n_state] Sigma_vec1;
  vector<lower = 30>[n_state] Sigma_vec2;
  vector<lower = 30>[n_state] Sigma_vec3;
  
  real<lower = 0> Sigma_tot;
  
  real<lower = 0, upper = 2> kk;
}

transformed parameters{
  
  matrix[n_state, n_state] alpha1;
  matrix[n_state, n_state] alpha2;
  matrix[n_state, n_state] alpha3;
  
  matrix[n_trans*n_state, n_trans*n_state] new_alpha1;
  matrix[n_trans*n_state, n_trans*n_state] new_alpha2;
  matrix[n_trans*n_state, n_trans*n_state] new_alpha3;
  
  alpha1 = lambda_to_alpha(knn_weight, n_state, lambda1, lambda_death);
  alpha2 = lambda_to_alpha(knn_weight, n_state, lambda2, lambda_death);
  alpha3 = lambda_to_alpha(knn_weight, n_state, lambda3, lambda_death);
  
  new_alpha1 = lambda_to_alpha_5(knn_weight, n_state, alpha1, alpha2, alpha3, s_rate, lambda_death, 0, 0, s_phase_ratio);
  new_alpha2 = lambda_to_alpha_5(knn_weight, n_state, alpha1, alpha2, alpha3, s_rate, lambda_death, lambda_trans1, 0, s_phase_ratio);
  new_alpha3 = lambda_to_alpha_5(knn_weight, n_state, alpha1, alpha2, alpha3, s_rate, lambda_death, 0, lambda_trans2, s_phase_ratio);
  //   
  // new_alpha1 = lambda_to_alpha_3(knn_weight, n_state, alpha1, alpha2, alpha3, lambda_death, 0, 0);
  // new_alpha2 = lambda_to_alpha_3(knn_weight, n_state, alpha1, alpha2, alpha3, lambda_death, lambda_trans1, 0);
  // new_alpha3 = lambda_to_alpha_3(knn_weight, n_state, alpha1, alpha2, alpha3, lambda_death, 0, lambda_trans2);
  
}

model {

  vector[n_trans*n_state] x_t_1;
  vector[n_trans*n_state] x_t_1_long;
  real x_t_1_long_sum;
  vector[n_state] x_t_0;
  
  vector[n_trans*n_state] x_0;
  vector[n_trans*n_state] x_11;
  vector[n_trans*n_state] x_15;
  
  // lambda_h ~ cauchy(0, 1);
  // tau ~ cauchy(0, 1);
  // for(i in 1:n_rate) lambda1[i] ~ normal(0, lambda_h[i] * tau);
  
  // lambda1     ~ normal(   0,    2);
  // lambda2     ~ normal(   0,    2);
  // lambda3     ~ normal(   0,    2);
  kk ~ normal(2,2);
  
  lambda1 ~ double_exponential(0, kk);
  lambda2 ~ double_exponential(0, kk);
  lambda3 ~ double_exponential(0, kk);
  
  lambda_death ~ normal(   0,    1);

  lambda_trans1  ~ normal(   0,    2);
  lambda_trans2  ~ normal(   0,    1);
  
  Sigma_vec1  ~ normal( 2000, 1500);
  Sigma_vec2  ~ normal( 2000, 1500);
  Sigma_vec3  ~ normal( 2000, 1500);
  
  Sigma_tot   ~ normal( 0, 10);
  
  s_rate ~ normal(0, 5);

  x_0 = rep_vector(0, n_trans*n_state);
  
  x_0[1:n_state] = x_t[1]*y_t_in[1];
  
  x_11  = matrix_exp(9 * new_alpha1) * x_0;
  x_15  = matrix_exp(4 * new_alpha2) * x_11;
  
  for(t in 2:n_time) {
    
      if(time_points[t]<11) x_t_1_long = matrix_exp((time_points[t]-time_points[1]) * new_alpha1) * x_0;
      else if(time_points[t]>=11 && time_points[t]<15) x_t_1_long = matrix_exp((time_points[t]-11) * new_alpha2) * x_11;
      else if(time_points[t]>=15) x_t_1_long = matrix_exp((time_points[t]-15) * new_alpha3) * x_15;
      
      // x_t_1_long = matrix_exp((time_points[t]-time_points[1]) * alpha) * x_0;
      x_t_1 = x_t_1_long[1:(n_trans*n_state)]+1e-6;
      x_t_1 = x_t_1/sum(x_t_1);
      
      for(m in 1:n_state){
      
      x_t_0[m] = x_t_1[m] + x_t_1[m+n_state] + x_t_1[m+2*n_state];
      
      }
        
      //   
      // x_t_0[1] = x_t_1[1] + x_t_1[5] + x_t_1[9];
      // x_t_0[2] = x_t_1[2] + x_t_1[6] + x_t_1[10];
      // x_t_0[3] = x_t_1[3] + x_t_1[7] + x_t_1[11];
      // x_t_0[4] = x_t_1[4] + x_t_1[8] + x_t_1[12];
      // 
      
      if(time_points[t]<11) x_t[t] ~ beta(Sigma_vec1 .* x_t_0, Sigma_vec1 .* (1-x_t_0));
      else if(time_points[t]>=11 && time_points[t]<15) x_t[t] ~ beta(Sigma_vec2 .* x_t_0, Sigma_vec2 .* (1-x_t_0));
      else if(time_points[t]>=15) x_t[t] ~ beta(Sigma_vec3 .* x_t_0, Sigma_vec3 .* (1-x_t_0));
    
      x_t_1_long_sum =  log10(sum(x_t_1_long));
     
      x_t_1_long_sum ~ normal(y_t[t], Sigma_tot);
      //print(sum(x_t_1_long));
   
  }
}

generated quantities{
  vector[n_trans*n_state] x_hat;
  vector[n_trans*n_state] x_0;
  vector[n_trans*n_state] x_11;
  vector[n_trans*n_state] x_15;
  
  vector[n_state] x_t_0;
  
  int t = 2;
  int total = 181;
  
  vector[n_state] x_t_1[181];
  vector[181] y_t_out;
  
  real x_sum;
  
  real new_time_point = 2.10;
  
  x_t_1[1] = to_vector(x_t[1]);
  
  y_t_out[1] = log10(y_t_in[1]);
    
    x_0 = rep_vector(0, n_trans*n_state);
    x_0[1:n_state] = x_t_1[1]*y_t_in[1];
    
    x_11 = matrix_exp((11-time_points[1]) * new_alpha1) * x_0;
    x_15 = matrix_exp(4 * new_alpha2) * x_11;
    
    while(t <= total){
      
      // print(new_time_point);
      
      if(new_time_point <11) x_hat = matrix_exp((new_time_point-time_points[1]) * new_alpha1) * x_0;
      if(new_time_point >=11 && new_time_point <15) x_hat = matrix_exp((new_time_point-11) * new_alpha2) * x_11;
      if(new_time_point >=15) x_hat = matrix_exp((new_time_point-15) * new_alpha3) * x_15;
      
      x_sum = sum(x_hat);
      
      x_hat = (x_hat + 1e-6)/sum(x_hat);
      
      for(m in 1:n_state){
        
        x_t_0[m] = x_hat[m] + x_hat[m+n_state] + x_hat[m+2*n_state];
        
      }

      //x_t_0[1] = x_hat[1] + x_hat[5] + x_hat[9];
      //x_t_0[2] = x_hat[2] + x_hat[6] + x_hat[10];
      //x_t_0[3] = x_hat[3] + x_hat[7] + x_hat[11];
      //x_t_0[4] = x_hat[4] + x_hat[8] + x_hat[12];
      
      // print(x_t_0);
      
      if(new_time_point <11) for(i in 1:n_state) x_t_1[t, i] = beta_rng(Sigma_vec1[i]* x_t_0[i], Sigma_vec1[i]* (1-x_t_0[i]));
      if(new_time_point >=11 && new_time_point <15) for(i in 1:n_state) x_t_1[t, i] = beta_rng(Sigma_vec2[i]* x_t_0[i], Sigma_vec2[i]* (1-x_t_0[i]));
      if(new_time_point >=15) {
        
        for(i in 1:n_state) {
          
          // if((Sigma_vec3[i]* (1-x_t_0[i]))<0) print(x_t_0);
          x_t_1[t, i] = beta_rng(Sigma_vec3[i]* x_t_0[i], Sigma_vec3[i]* (1-x_t_0[i]));
          
        }
      }
      
      //from this: sum(x_t_1_long) ~ normal(y_t[t], Sigma_tot);
      y_t_out[t] = normal_rng(log10(x_sum), Sigma_tot);
      //y_t_out[t] = x_sum;
      
      t = t + 1;
      new_time_point = new_time_point + 0.10;
    }
}

//  generated quantities{
//   vector[n_rate*n_state] x_hat;
//   vector[n_rate*n_state] x_0;
//   vector[n_rate*n_state] x_11;
//   vector[n_rate*n_state] x_15;
// 
//   vector[n_state] x_t_0;
// 
//   int t = 2;
//   int total = 181;
// 
//   // matrix[n_state, n_state] alpha1;
//   // matrix[n_state, n_state] alpha2;
//   // matrix[n_state, n_state] alpha3;
//   //
//   // matrix[n_rate*n_state, n_rate*n_state] new_alpha1;
//   // matrix[n_rate*n_state, n_rate*n_state] new_alpha2;
//   // matrix[n_rate*n_state, n_rate*n_state] new_alpha3;
// 
//   vector[n_state] x_t_1[181];
//   vector[181] y_t_out;
// 
//   real x_sum;
// 
//   real new_time_point = 2.10;
// 
//   x_t_1[1] = to_vector(x_t[1]);
// 
//   y_t_out[1] = log10(y_t_in[1]);
// 
//   // alpha1 = lambda_to_alpha_3(knn_weight, n_state, lambda1, lambda_death);
//   // alpha2 = lambda_to_alpha_3(knn_weight, n_state, lambda2, lambda_death);
//   // alpha3 = lambda_to_alpha_3(knn_weight, n_state, lambda3, lambda_death);
//   //
//   // new_alpha1 = lambda_to_alpha_5(knn_weight, n_state, alpha1, alpha2, alpha3, s_rate, lambda_death, 0, 0, s_phase_ratio);
//   // new_alpha2 = lambda_to_alpha_5(knn_weight, n_state, alpha1, alpha2, alpha3, s_rate, lambda_death, lambda_trans1, 0, s_phase_ratio);
//   // new_alpha3 = lambda_to_alpha_5(knn_weight, n_state, alpha1, alpha2, alpha3, s_rate, lambda_death, 0, lambda_trans2, s_phase_ratio);
// 
//   x_0 = rep_vector(0, n_rate*n_state);
//   x_0[1:4] = x_t_1[1]*y_t_in[1];
// 
//   x_11 = matrix_exp((11-time_points[1]) * new_alpha1) * x_0;
//   x_15 = matrix_exp(4 * new_alpha2) * x_11;
// 
//   // print(x_11);
//   // print(x_15);
// 
//   while(t <= total){
// 
//     // print(new_time_point);
// 
//     if(new_time_point <11) x_hat = matrix_exp((new_time_point-time_points[1]) * new_alpha1) * x_0;
//     if(new_time_point >=11 && new_time_point <15) x_hat = matrix_exp((new_time_point-11) * new_alpha2) * x_11;
//     if(new_time_point >=15) x_hat = matrix_exp((new_time_point-15) * new_alpha3) * x_15;
// 
//     x_sum = sum(x_hat);
// 
//     x_hat = (x_hat + 1e-6)/sum(x_hat);
// 
//     x_t_0[1] = x_hat[1] + x_hat[5] + x_hat[9];
//     x_t_0[2] = x_hat[2] + x_hat[6] + x_hat[10];
//     x_t_0[3] = x_hat[3] + x_hat[7] + x_hat[11];
//     x_t_0[4] = x_hat[4] + x_hat[8] + x_hat[12];
// 
//     // print(x_t_0);
// 
//     if(new_time_point <11) for(i in 1:n_state) x_t_1[t, i] = beta_rng(Sigma_vec1[i]* x_t_0[i], Sigma_vec1[i]* (1-x_t_0[i]));
//     if(new_time_point >=11 && new_time_point <15) for(i in 1:n_state) x_t_1[t, i] = beta_rng(Sigma_vec2[i]* x_t_0[i], Sigma_vec2[i]* (1-x_t_0[i]));
//     if(new_time_point >=15) {
// 
//       for(i in 1:n_state) {
// 
//         // if((Sigma_vec3[i]* (1-x_t_0[i]))<0) print(x_t_0);
//         x_t_1[t, i] = beta_rng(Sigma_vec3[i]* x_t_0[i], Sigma_vec3[i]* (1-x_t_0[i]));
// 
//     }
//     }
// 
//     //from this: sum(x_t_1_long) ~ normal(y_t[t], Sigma_tot);
//     y_t_out[t] = normal_rng(log10(x_sum), Sigma_tot);
//     //y_t_out[t] = x_sum;
// 
//     t = t + 1;
//     new_time_point = new_time_point + 0.10;
//   }
// }
