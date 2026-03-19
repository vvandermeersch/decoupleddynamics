functions {
  vector gp_pred_rng(array[] real x2,
                     vector y1, array[] real x1,
                     vector mu,
                     real alpha, real rho,
                     real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   gp_exp_quad_cov(x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1 - mu);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = gp_exp_quad_cov(x1, x2, alpha, rho);
      vector[N2] f2_mu = mu + (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   gp_exp_quad_cov(x2, alpha, rho)
                              - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  int<lower=1> N; 
  int<lower=1> N_weeks; 
  array[N] int days; 
  array[N_weeks] int weeks;
  array[N] int week_idxs; 
  vector[N] nee_obs;
}

parameters {
  real alpha;         
  real<lower=0> rho;   
  real<lower=0> gamma; 
  real<lower=0> sigma; 
  
  vector[N_weeks] fweek_tilde;  
}

transformed parameters{
  matrix[N_weeks, N_weeks] cov =  gp_exp_quad_cov(weeks, gamma, rho)
                       + diag_matrix(rep_vector(1e-10, N_weeks));
  matrix[N_weeks, N_weeks] L_cov = cholesky_decompose(cov);
  
  vector[N_weeks] fweek = alpha + L_cov * fweek_tilde; 
}

model {
  alpha ~ normal(0, 1);
  rho ~ normal(50, 10);  
  gamma ~ normal(0, 1); 
  sigma ~ normal(0, 1);

  fweek_tilde ~ normal(0,1);
  // fweek ~ multi_normal_cholesky(rep_vector(alpha, N_weeks), L_cov);
  
  nee_obs ~ normal(fweek[week_idxs], sigma);
}

generated quantities {
  
  array[N] real nee_pred;
  
  nee_pred = normal_rng(fweek[week_idxs], sigma);

}
