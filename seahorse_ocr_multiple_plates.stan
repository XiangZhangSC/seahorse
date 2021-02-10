data {
  int<lower=0> N_obs;
  vector<lower=0>[N_obs] OCR_obs;
  int<lower=0> N_lev1; // Number of combinations Cell, Injection, Plate and Well
  int<lower=1,upper=N_lev1> id1[N_obs];
  int<lower=0> N_lev2; // Number of combinations Cell, Injection and Plate
  int<lower=1,upper=N_lev2> id2[N_lev1];
  int<lower=0> N_lev3; // Number of combinations Cell and Injection
  int<lower=1,upper=N_lev3> id3[N_lev2];
  int<lower=0> N_lev4; // Number of combinations Plate and Injection
  int<lower=1,upper=N_lev4> id4[N_obs];
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
  positive_ordered[4] control_OCR_per_1k;
  positive_ordered[4] patient_OCR_per_1k;
  vector<lower=0>[N_lev2] OCR_per_1k;
  vector[N_lev1] mulog_OCR; // the magnitude of OCR in a particular well and with a particular injection
}
transformed parameters {
  vector<lower=0>[N_lev3] mu_OCR_per_1k;
  
  mu_OCR_per_1k[1] = control_OCR_per_1k[3]; // initial
  mu_OCR_per_1k[2] = control_OCR_per_1k[2]; // oligomycin
  mu_OCR_per_1k[3] = control_OCR_per_1k[4]; // fccp
  mu_OCR_per_1k[4] = control_OCR_per_1k[1]; // rotenone
  
  mu_OCR_per_1k[5] = patient_OCR_per_1k[3]; // initial
  mu_OCR_per_1k[6] = patient_OCR_per_1k[2]; // oligomycin
  mu_OCR_per_1k[7] = patient_OCR_per_1k[4]; // fccp
  mu_OCR_per_1k[8] = patient_OCR_per_1k[1]; // rotenone
}
model {
  
  sdlog_OCR ~ lognormal( -3.23, 0.79 );
  sigma_well ~ lognormal( -1.62, 0.53 );
  sigma_plate ~ lognormal( -1.18, 0.05 );
  
  // prior for initial phase
  control_OCR_per_1k[3] ~ lognormal( 0.75, 0.29 );
  patient_OCR_per_1k[3] ~ lognormal( 0.75, 0.29 );
  
  // prior for oligomycin phase
  control_OCR_per_1k[2] ~ lognormal( -0.28, 0.32 );
  patient_OCR_per_1k[2] ~ lognormal( -0.28, 0.32 );
  
  // prior for FCCP phase
  control_OCR_per_1k[4] ~ lognormal( 1.26, 0.30 );
  patient_OCR_per_1k[4] ~ lognormal( 1.26, 0.30 );
  
  
  // prior for rotenone phase
  control_OCR_per_1k[1] ~ lognormal( -0.74, 0.32 );
  patient_OCR_per_1k[1] ~ lognormal( -0.74, 0.32 );
  
  for (k in 1:N_lev2 ) {
    OCR_per_1k[k] ~ lognormal(log(mu_OCR_per_1k[id3[k]]), sigma_plate[id3[k]]);
  }
  
  for ( j in 1:N_lev1 ) {
    mulog_OCR[j] ~ normal(log(OCR_per_1k[id2[j]] * N_cell_1k[j]), sigma_well[id2[j]]);
  }
  
  // likelihood 
  for ( i in 1:N_obs ) {
    OCR_obs[i] ~ lognormal( mulog_OCR[id1[i]], sdlog_OCR[id4[i]] );
  }
}
