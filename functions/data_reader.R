#  data_reader.R
#  Created by John T. Trochta
#  This file reads in .DAT and .CTL files created for ADMB .TPL files
#  The use must ensure not to have within line comments in the files (nope, fixed that)

data.reader <- function(filename) {
    #  The user needs to make sure there is a blank line at the end of
    #  filename and that each data type (vector number or matrix) is
    #  separated by a blank line
    
    # This is kind of convoluted
    text <- readLines(filename)
    values <- grep("^\\s{0,2}[0-9]", text)
    signed.values <- grep("^\\s{0,2}[-]", text)
    read.these <- sort(c(values, signed.values))
    nlines <- length(text) 
    indices <- seq(1:nlines)
    indices <- indices[read.these]
    first.differences <- c(diff(indices),5)# This accounts for the last data
    data.types <- length(first.differences[first.differences>1])
    
    data <- vector("list", data.types)
    j <- 1
    temp <- NA
    for(i in 1:length(indices)){
        temp.1 <- scan(filename, skip=indices[i]-1, nlines=1, quiet=TRUE, flush=FALSE)
        if(first.differences[i]>1 | all(first.differences[1:(length(first.differences)-1)]==1)){
            if(all(is.na(temp))){
                data[[j]] <- temp.1
                j <- j+1
            }else{
                temp <- rbind(temp, temp.1)
                data[[j]] <- temp
                rownames(data[[j]]) <- NULL
                temp <- NA
                j <- j+1
            }
        } else{
            if(all(is.na(temp))){
                temp <- temp.1
            }else{
                temp <- rbind(temp, temp.1)
            }
        }
    }
    return(data)
}

# Read BASA data files (.dat and .ctl) into list format for easy parameter access.
read.data.files <- function(dat.dir){

    filename <- vector(length=5)
    filename[1] <- paste0(dat.dir, "PWS_ASA.dat")
    filename[2] <- paste0(dat.dir, "PWS_ASA(ESS).ctl")
    filename[3] <- paste0(dat.dir, "PWS_ASA(covariate).ctl")
    filename[4] <- paste0(dat.dir, "agecomp_samp_sizes.txt")
    filename[5] <- paste0(dat.dir, "PWS_ASA_disease.dat")

    PWS_ASA.dat            <- data.reader(filename=filename[1])
    PWS_ASA_ESS.ctl        <- data.reader(filename=filename[2])
    PWS_ASA_covariate.ctl  <- data.reader(filename=filename[3])
    agecomp_samp_sizes.txt <- data.reader(filename=filename[4])
    PWS_ASA_disease.dat    <- data.reader(filename=filename[5])

    par.names <- c("nyr", "nyr_tobefit", "nage", "waa", "fecundity", 
                "pound_catch", "pk", "foodbait_catch", "gillnet_catch", "seine_yield",  "perc.female", 
                "mdm", "egg", "egg_se", "adfg_hydro_year_start", "adfg_hydro", "pwssc_hydro_year_start", "pwssc_hydro", "pwssc_hydro_se",
                "seine_age_comp", "spawn_age_comp", "juvenile_survey")
    names(PWS_ASA.dat) <- par.names

    par.names <- c("ctl_file", "seine_ess", "spawn_ess", "vhsv_ess", "ich_ess")
    names(PWS_ASA_ESS.ctl) <- par.names

    par.names <- c("std_covs", "num_recruit_cov", "recruit_fixed_rand", "cov_on", "regime_shift_89", "r_beta_change", 
                    "num_mort_covs", "mort_fixed_rand", "mort_on", "mort_age_impact", "disease_covs", "winter_mort_devs", "m_btea_change")
    names(PWS_ASA_covariate.ctl) <- par.names

    par.names <- c("seine_sample_size", "spawn_sample_size", "vhsv_sample_size", "ich_smaple_size")
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

# Function for simultaneously creating names of variables/elements within a list
listN <- function(...){
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    return(anonList)
}
