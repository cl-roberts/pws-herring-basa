#' Read Data Files
#'
#' Read in the five data files containing survey data.
#'
#' @param dat.dir A character string giving the directory containg the `.dat`,
#' `.ctl`, and `.txt` data files
#'
#' @returns A named list with an element corresponding to each data file.


read.data.files <- function(dat.dir){

  # CLR: add / to file paths
  filename <- vector(length=5)
  filename[1] <- paste0(dat.dir, "/PWS_ASA.dat")
  filename[2] <- paste0(dat.dir, "/PWS_ASA(ESS).ctl")
  filename[3] <- paste0(dat.dir, "/PWS_ASA(covariate).ctl")
  filename[4] <- paste0(dat.dir, "/agecomp_samp_sizes.txt")
  filename[5] <- paste0(dat.dir, "/PWS_ASA_disease.dat")

  PWS_ASA.dat            <- data.reader(filename=filename[1])
  PWS_ASA_ESS.ctl        <- data.reader(filename=filename[2])
  PWS_ASA_covariate.ctl  <- data.reader(filename=filename[3])
  agecomp_samp_sizes.txt <- data.reader(filename=filename[4])
  PWS_ASA_disease.dat    <- data.reader(filename=filename[5])

  par.names <- c("nyr", "nyr_tobefit", "nage", "waa", "fecundity",
                 "pound_catch", "pk", "foodbait_catch", "gillnet_catch", "seine_yield",  "perc_female",
                 "mdm", "egg", "egg_se", "adfg_hydro_year_start", "adfg_hydro", "pwssc_hydro_year_start", "pwssc_hydro", "pwssc_hydro_se",
                 "seine_age_comp", "spawn_age_comp", "juvenile_survey")
  names(PWS_ASA.dat) <- par.names

  par.names <- c("ctl_file", "seine_ess", "spawn_ess", "vhsv_ess", "ich_ess")
  names(PWS_ASA_ESS.ctl) <- par.names

  par.names <- c("std_covs", "num_recruit_cov", "recruit_fixed_rand", "cov_on", "regime_shift_89", "r_beta_change",
                 "num_mort_covs", "mort_fixed_rand", "mort_season", "mort_on", "mort_age_impact", "disease_covs", "winter_mort_devs", "m_beta_change")
  names(PWS_ASA_covariate.ctl) <- par.names

  par.names <- c("seine_sample_size", "spawn_sample_size", "vhsv_sample_size", "ich_sample_size")
  names(agecomp_samp_sizes.txt) <- par.names

  par.names <- c("vhsv_age_prevalence", "vhsv_obs_start", "vhsv_est_start", "vhsv_recov_prob",
                 "ich_age_prevalence", "ich_obs_start", "ich_est_start", "ich_recov_prob")
  names(PWS_ASA_disease.dat) <- par.names

  return(listN(PWS_ASA.dat,
               PWS_ASA_ESS.ctl,
               PWS_ASA_covariate.ctl,
               agecomp_samp_sizes.txt,
               PWS_ASA_disease.dat))
}
