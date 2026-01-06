update_om_data <- function(data, ghl, new_waa, new_perc_female, new_disease_covs, new_seine_ess, new_spawn_ess, KI_ghl = NULL, new_KI_seine_ess = NULL, new_KI_spawn_ess = NULL) {

    om_data <- data

    neg_9s <- rep(-9, om_data$nage)
    zeros  <- rep(0, om_data$nage)

    # PWS_ASA.dat
    om_data$nyr <- om_data$nyr+1
    om_data$waa <- rbind(om_data$waa, new_waa)
    om_data$fecundity <- rbind(om_data$fecundity, neg_9s)
    om_data$pound_catch <- rbind(om_data$pound_catch, zeros)
    om_data$foodbait_catch <- rbind(om_data$foodbait_catch, zeros)
    om_data$gillnet_catch <- rbind(om_data$gillnet_catch, zeros)
    om_data$seine_yield <- rbind(om_data$seine_yield, ghl)
    om_data$perc_female <- rbind(om_data$perc_female, new_perc_female)
    om_data$mdm <- rbind(om_data$mdm, -9)
    om_data$egg <- rbind(om_data$egg, -9)
    om_data$egg_se <- rbind(om_data$egg_se, -9)
    om_data$adfg_hydro <- rbind(om_data$adfg_hydro, -9)
    om_data$pwssc_hydro <- rbind(om_data$pwssc_hydro, -9)
    om_data$pwssc_hydro_se <- rbind(om_data$pwssc_hydro_se, -9)
    om_data$seine_age_comp <- rbind(om_data$seine_age_comp, neg_9s)
    om_data$spawn_age_comp <- rbind(om_data$spawn_age_comp, neg_9s)
    om_data$juvenile_survey <- rbind(om_data$juvenile_survey, -9)

    # PWS_ASA_covariate
    om_data$disease_covs <- rbind(om_data$disease_covs, new_disease_covs)

    # agecomp sample sizes
    om_data$seine_sample_size <- rbind(om_data$seine_sample_size, 0)
    om_data$spawn_sample_size <- rbind(om_data$spawn_sample_size, 0)

    # PWS_ASA_ESS
    om_data$seine_ess <- rbind(om_data$seine_ess, new_seine_ess)
    om_data$spawn_ess <- rbind(om_data$spawn_ess, new_spawn_ess)

    if (!is.null(KI_ghl)) {
        om_data$KI_seine_yield <- rbind(om_data$KI_seine_yield, KI_ghl)
    }
    
    if (!is.null(new_KI_seine_ess)) {
        om_data$KI_seine_ess <- rbind(om_data$KI_seine_ess, new_KI_seine_ess)
    }
    
    if (!is.null(new_KI_seine_ess)) {
        om_data$KI_seine_age_comp <- rbind(om_data$KI_seine_age_comp, neg_9s)
    }

    if (!is.null(new_KI_spawn_ess)) {
        om_data$KI_spawn_ess <- rbind(om_data$KI_spawn_ess, new_KI_spawn_ess)
    }
    

    return(om_data)

}