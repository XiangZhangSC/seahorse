library(dplyr)
library(tidybayes)
library(broom)
library(ggplot2)
library(forcats)

find_host_plate <- function(df, target_cell_line) {
  df %>% 
    dplyr::select(plate_id, Cell_line) %>% 
    distinct() %>% 
    dplyr::filter(Cell_line %in% target_cell_line) %>% 
    dplyr::pull(plate_id)
}

aggregate_measurements <- function(df) {
  # three measurement cycles within an interval
  dat.lev1 <- df %>% 
    dplyr::select(Cell_line, Interval, plate_id, well, cell_n) %>% 
    distinct()
  
  dat.lev1$id1 <- row.names(dat.lev1)
  
  return(dat.lev1)
}

aggregate_wells <- function(df, single_plate = TRUE) {
  
  dat.lev2 <- df %>% 
      dplyr::select(Cell_line, Interval, plate_id) %>% 
      distinct()
  
  dat.lev2$id2 <- row.names(dat.lev2)
  
  return(dat.lev2)
}

aggregate_plates <- function(df) {
  
  dat.lev3 <- df %>% 
      dplyr::select(Cell_line, Interval) %>% 
      distinct()
    
  dat.lev3$id3 <- row.names(dat.lev3)
  
  return(dat.lev3)
}

plate_interval <- function(df) {
  dat.lev4 <- df %>% 
    dplyr::select(plate_id, Interval) %>% 
    distinct()
  
  dat.lev4$id4 <- row.names(dat.lev4)
  
  return(dat.lev4)
}

collect_data_for_stan <- function(df) {
  df <- sort_case_dat(df)
  dat.lev1 <- aggregate_measurements(df)
  dat.lev2 <- aggregate_wells(dat.lev1)
  dat.lev3 <- aggregate_plates(dat.lev2)
  dat.lev4 <- plate_interval(df)
  
  dat.lev0.new <- df %>% 
    left_join(dat.lev1) %>% 
    left_join(dat.lev4)
  
  dat.lev1.new <- dat.lev1 %>% 
    dplyr::select(-id1) %>% 
    left_join(dat.lev2)
  
  dat.lev2.new <- dat.lev2 %>% 
    dplyr::select(-id2) %>% 
    left_join(dat.lev3)
  
  dat_list_for_stan <- list(
    N_obs = nrow(dat.lev0.new), 
    OCR_obs = dat.lev0.new$OCR, 
    N_lev1 = nrow(dat.lev1.new), 
    id1 = as.integer(dat.lev0.new$id1), 
    N_lev2 = nrow(dat.lev2.new), 
    id2 = as.integer(dat.lev1.new$id2), 
    N_lev3 = nrow(dat.lev3), 
    id3 = as.integer(dat.lev2.new$id3),
    N_lev4 = nrow(dat.lev4),
    id4 = as.integer(dat.lev0.new$id4),
    N_cell = dat.lev1.new$cell_n
  )
  
  return(dat_list_for_stan)
}

tabulate_posterior_ocr <- function(stanfit.obj, experiment_setup) {
  post_ocr_1k <- spread_draws(stanfit.obj, mu_OCR_per_1k[id3])
  
  post_ocr_1k_final <- post_ocr_1k %>% 
    mutate(id3 = as.character(id3)) %>% 
    left_join(experiment_setup, by = "id3")
  
  return(post_ocr_1k_final)
}

convert_ocr_to_respiration <- function(stanfit.obj, experiment_setup) {
  
  post_ocr_1k_final <- tabulate_posterior_ocr(stanfit.obj, experiment_setup)
  
  post_respiration.df <- post_ocr_1k_final %>% 
    ungroup() %>% 
    dplyr::select(-id3) %>% 
    spread(Interval, mu_OCR_per_1k) %>% 
    mutate(`Basal respiration` = Int1 - Int4, 
           `ATP-linked respiration` = Int1 - Int2, 
           `Proton leak` = Int2 - Int4, 
           `Spare respiratory capacity` = Int3 - Int1, 
           `Maximal respiration` = Int3 - Int4, 
           Int1 = NULL, 
           Int2 = NULL, 
           Int3 = NULL, 
           Int4 = NULL) %>% 
    gather(what, val, `Basal respiration`:`Maximal respiration`)
}

sketch_respirations <- function(df, my_title) {
  df %>% 
    mutate(what == fct_inorder(what)) %>% 
    ggplot(aes(val)) + 
    geom_density(aes(color = Cell_line)) + 
    scale_color_brewer("Cell line", palette = "Set2") + 
    facet_wrap(~what, ncol = 1, scales = "free_y") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    labs(x = "Respiration", title = my_title)
}

true_signal_measurement_error <- function(x) {
  tidy(MASS::fitdistr(x, densfun = "lognormal"))
}

collect_data_per_case <- function(patient_dat, control_dat) {
  
  patient_dat$Cell_line <- "patient"
  
  related_plates <- unique(patient_dat$plate_id)
  control_dat <- control_dat %>% 
    dplyr::filter(plate_id %in% related_plates)
  
  case_dat <- control_dat %>% 
    bind_rows(patient_dat)
  
  return(case_dat)
}

sort_case_dat <- function(case_dat) {
  case_dat$Cell_line <- as.factor(case_dat$Cell_line)
  case_dat$Cell_line <- relevel(case_dat$Cell_line, ref = "NHDF")
  case_dat_sorted <- case_dat %>% 
    arrange(plate_id, Cell_line, Interval)
  return(case_dat_sorted)
}

tag_experimental_setup <- function(case_dat) {
  case_dat$Cell_line <- as.factor(case_dat$Cell_line)
  case_dat$Cell_line <- relevel(case_dat$Cell_line, ref = "NHDF")
  case_dat <- case_dat %>% 
    arrange(Cell_line)
  dat.lev1 <- aggregate_measurements(case_dat)
  dat.lev2 <- aggregate_wells(dat.lev1)
  dat.lev3 <- aggregate_plates(dat.lev2)
  return(dat.lev3)
}

experimental_setup_per_case <- function(patient_dat, control_dat) {
  
  case_dat <- collect_data_per_case(patient_dat, control_dat)
  experimental_setup <- tag_experimental_setup(case_dat)
  
  return(experimental_setup)
} 

plate_phase_per_case <- function(patient_dat, control_dat) {
  
  case_dat <- collect_data_per_case(patient_dat, control_dat)
  dat.lev4 <- plate_interval(case_dat)
  
  return(dat.lev4)
} 

run_ocrbayes_per_case <- function(patient_dat, control_dat) {
  case_dat <- collect_data_per_case(patient_dat, control_dat)
  case_fit <- ocrbayes(case_dat)
  return(case_fit)
}

add_fdr <- function(df) {
  my_fdrs <- vector("numeric", length = nrow(df))
  df.sorted <- df %>% arrange(PEP)
  for (i in 1:nrow(df)) {
    my_fdrs[i] <- mean(head(df.sorted$PEP, i))
  }
  return(my_fdrs)
}

add_var_cycles <- function(patient_dat, control_dat, stanobj) {
  my_var_cycles <- spread_draws(stanobj, sdlog_OCR[id4])
  
  my_plate_phase <- plate_phase_per_case(patient_dat, control_dat)
  
  my_var_cycles_final <- my_var_cycles %>% 
    mutate(id4 = as.character(id4)) %>% 
    left_join(my_plate_phase)
  
  return(my_var_cycles_final)
}
