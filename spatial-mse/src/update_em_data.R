update_em_data <- function(data, ghl, new_waa, new_perc_female, new_disease_covs, new_seine_ess, new_spawn_ess, new_mdm, new_juv, new_seine_agecomp, new_spawn_agecomp, new_seine_sample_size, new_spawn_sample_size) {

    em_data <- update_om_data(
        data, ghl = ghl, new_waa = new_waa, new_perc_female = new_perc_female,
        new_disease_covs = new_disease_covs, new_seine_ess = new_seine_ess, new_spawn_ess = new_spawn_ess, 
        KI_ghl = NULL, new_KI_seine_ess = NULL, new_KI_spawn_ess = NULL
    )

    # PWS_ASA.dat
    em_data$mdm[em_data$nyr] <- new_mdm
    em_data$seine_age_comp[em_data$nyr,] <- new_seine_agecomp
    em_data$spawn_age_comp[em_data$nyr,] <- new_spawn_agecomp
    em_data$juvenile_survey[em_data$nyr] <- new_juv

    # raw sample sizes
    em_data$seine_sample_size[em_data$nyr] <- new_seine_sample_size  
    em_data$spawn_sample_size[em_data$nyr,] <- new_spawn_sample_size

    return(em_data)

}