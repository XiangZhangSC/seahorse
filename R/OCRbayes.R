source("./R/customized_functions.R")
library(rstan)

my_stan_mod <- stan_model("seahorse_ocr_multiple_plates.stan")

ocrbayes <- function(df, iter = 2000, num.chains = 4L, num.threads = 4L, max_treedepth = 10, adapt_delta = 0.95) {
  dat_list_for_stan <- collect_data_for_stan(df)
  stan_output <- sampling(my_stan_mod, data = dat_list_for_stan, 
                          iter = iter, 
                          chains = num.chains, cores = num.threads, 
                          control = list(max_treedepth = max_treedepth, adapt_delta = adapt_delta))
  return(stan_output)
}

