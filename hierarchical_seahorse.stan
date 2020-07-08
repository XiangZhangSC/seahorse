//
// This Stan program defines a hierachical model for a
// 96-well Seahorse XF plate.
//

data {
  int<lower=0> N;
  vector[N] y_obs;
  int<lower=0> N_well;
  int<lower=1, upper=N_well> well[N];
  int<lower=0> N_injection;
  int<lower=1, upper=N_injection> injection[N];
}

parameters {
  real<lower=0> sigma;
  real z[N_well];
  real b[N_injection];
  real<lower=0> a_bar;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y[N_injection];
  matrix[N_well, N_injection] g;
}

model {
  matrix[N_well, N_injection] y_true;
  
  a_bar ~ normal( 0, 100 );
  sigma_a ~ exponential( 1 );
  z ~ normal( 0, 1 );
  to_vector(g) ~ normal( 0, 1 );
  b ~ normal( 0, 2.5 );
  sigma ~ exponential( 1 );
  sigma_y ~ exponential( 1 );
  
  for ( i in 1:N ) {
    // y_true[i] ~ normal( mu[i], sigma )
    // Biological OCR + injection effect
    y_true[well[i], injection[i]] = a_bar + z[well[i]] * sigma_a + b[injection[i]] + g[well[i], injection[i]] * sigma;
    
    // Every injection has three measurements
    y_obs[i] ~ normal(y_true[well[i], injection[i]], sigma_y[injection[i]] );
  }
}
generated quantities {
  vector[N_well] basal_ocr;
  matrix[3, N_injection] y_pred;
  
  for ( j in 1:N_well ) {
    basal_ocr[j] = a_bar + z[j] * sigma_a;
  }
  
  for ( m in 1:3 ) {
    for ( n in 1:N_injection ) {
      y_pred[m, n] = normal_rng(basal_ocr[1] + b[n], sigma_y[n]);
    }
  }
}
