#' Run BASA on peeled data
#'
#' Helper function to run BASA on peeled data for retrospective analysis. Peeled
#' data is returned by \link[pwsHerringBasa]{write_retro_data}.
#'
#' @param i Number of years to peel
#' @param dir_model Character string giving model directory
#' @param dir_retro Character string giving retrospective directory
#'
#' @returns A copy of (peeled) BASA outputs in `dir_retro`.
#'

run.basa.retro <- function(i, dir_model, dir_retro){
    
    dir_retro_i <- here::here(dir_retro, paste0("basa-", i))
    if(!dir.exists(dir_retro_i)){
        dir.create(dir_retro_i)
    }

    if(!file.exists(here::here(dir_retro_i, "mcmc_out/PFRBiomass.csv"))){

        print(paste("mcmc_out directory does not exist in", dir_retro_i))

        file.copy(here::here(dir_model, "PWS_ASA.TPL"), dir_retro_i)
        file.copy(here::here(dir_model, "PWS_ASA(par).ctl"), dir_retro_i)
        file.copy(here::here(dir_model, "PWS_ASA.dat"), dir_retro_i)
        file.copy(here::here(dir_model, "agecomp_samp_sizes.txt"), dir_retro_i)    
        file.copy(here::here(dir_model, "PWS_ASA(covariate).ctl"), dir_retro_i)            
        file.copy(here::here(dir_model, "PWS_ASA_disease.dat"), dir_retro_i)

        dat.files <- read.data.files(dir_model)
        write_retro_dat(dat.files, i, dir_retro_i)
        run.basa(dir_retro_i)    

    } else {

        print("mcmc_out directory already exists")

    }
}
