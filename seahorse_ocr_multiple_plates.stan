data {
  int<lower=0> N_obs;
  vector<lower=0>[N_obs] OCR_obs;
  int<lower=0> N_lev1; // Number of combinations Cell, Injection, Plate and Well
  int<lower=1,upper=N_lev1> id1[N_obs];
  int<lower=0> N_lev2; // Number of combinations Cell, Injection and Plate
  int<lower=1,upper=N_lev2> id2[N_lev1];
  int<lower=0> N_lev3; // Number of combinations Cell and Injection
  int<lower=1,upper=N_lev3> id3[N_lev2];
  int<lower=0> N_lev4;
  int<lower=1,upper=N_lev4> id4[N_obs]; // combinations Plate and Injection
  vector<lower=0>[N_lev1] N_cell;
}
transformed data {
  vector<lower=0>[N_lev1] N_cell_1k;
  
  N_cell_1k = N_cell / 1000;
}
parameters {
  vector<lower=0>[N_lev4] sdlog_OCR; // variation in the magnitude of OCR between the three measurement cycles
  vector<lower=0>[N_lev2] sigma_well; // variability between replicate wells after accounting cell number difference
  vector<lower=0>[N_lev3] sigma_plate; // variability between replicate plates
  vector<lower=0>[N_lev3] mu_OCR_per_1k;
  vector[N_lev1] z1;
  vector[N_lev2] z2;
}
transformed parameters {
  vector<lower=0>[N_lev2] OCR_per_1k;
  vector[N_lev1] mulog_OCR; // the magnitude of OCR in a particular well and with a particular injection
  
  for (k in 1:N_lev2 ) {
    OCR_per_1k[k] = exp(log(mu_OCR_per_1k[id3[k]]) + sigma_plate[id3[k]] * z2[k]);
  }
  
  for ( j in 1:N_lev1 ) {
    mulog_OCR[j] = log(OCR_per_1k[id2[j]] * N_cell_1k[j]) + sigma_well[id2[j]] * z1[j];
  }
}
model {
  sdlog_OCR ~ lognormal( -3.23, 0.79 );
  sigma_well ~ lognormal( -1.62, 0.53 );
  sigma_plate ~ lognormal( -1.18, 0.05 );
  mu_OCR_per_1k ~ lognormal( 0.3, 0.79 );
  z1 ~ normal( 0, 1 );
  z2 ~ normal( 0, 1 );
  
  for ( i in 1:N_obs ) {
    OCR_obs[i] ~ lognormal( mulog_OCR[id1[i]], sdlog_OCR[id4[i]] );
  }
}
