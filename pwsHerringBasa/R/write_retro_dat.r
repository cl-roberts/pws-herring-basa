#' Write peeled data for retrospective analysis
#'
#' Peels terminal years off of BASA input data time series and writes to `.dat`
#' files. This peeled data is then fed to \link[pwsHerringBasa]{run.basa.retro}
#' which copies it to the correct directory and runs BASA. 
#'
#' @param dat.files an R object containing BASA data files (returned by 
#' \link[pwsHerringBasa]{read.data.files})
#' @param years.to.remove Integer giving the number of years to peel
#' @param dir_retro Character string giving retrospective directory
#'
#' @returns A copy of (peeled) BASA inputs in `dir_retro`.
#'

write_retro_dat <- function(dat.files, years.to.remove, dir_retro){

    write.data <- function(data, fname){
        write.table(
            rbind(" ", data[[1]], " "),
            file = fname, 
            append = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE
        )
        for(i in 2:length(data)){
            write.table(
                rbind(data[[i]], " "),
                file = fname, 
                append = TRUE, sep = " ",
                row.names = FALSE, col.names = FALSE,
                quote = FALSE
            )
        }
    }

    options(scipen=10)

    # Write main data file with actual catches (known w/ error)
    # and observed survey values (from fun_obsm.R).
    PWS.ASA.dat <- dat.files$PWS_ASA.dat
    # Add new data 
    PWS.ASA.dat[[1]]  <- PWS.ASA.dat[[1]]-years.to.remove                                              
    # Number of year in data
    PWS.ASA.dat[[2]]  <- PWS.ASA.dat[[2]]-years.to.remove                                              
    # Number of years to be fit to
    # 3 is number of age classes (10). Needs no change.
    PWS.ASA.dat[[4]]  <- PWS.ASA.dat[[4]][1:(nrow(PWS.ASA.dat[[4]])-years.to.remove),]                 
    # Weight-at-age 
    PWS.ASA.dat[[5]]  <- PWS.ASA.dat[[5]][1:(nrow(PWS.ASA.dat[[5]])-years.to.remove),]                 
    # Fecundity
    PWS.ASA.dat[[6]]  <- PWS.ASA.dat[[6]][1:(nrow(PWS.ASA.dat[[6]])-years.to.remove),]                 
    # Pound catch fishery
    # 7 is prop of pound-catch killed (0.25). Needs no change.
    PWS.ASA.dat[[8]]  <- PWS.ASA.dat[[8]][1:(nrow(PWS.ASA.dat[[8]])-years.to.remove),]                 
    # Foodbait fishery
    PWS.ASA.dat[[9]]  <- PWS.ASA.dat[[9]][1:(nrow(PWS.ASA.dat[[9]])-years.to.remove),]                 
    # Gillnet fishery
    PWS.ASA.dat[[10]] <- as.matrix(PWS.ASA.dat[[10]][1:(nrow(PWS.ASA.dat[[10]])-years.to.remove),])    
    # Seine net fishery
    PWS.ASA.dat[[11]] <- as.matrix(PWS.ASA.dat[[11]][1:(nrow(PWS.ASA.dat[[11]])-years.to.remove),])    
    # % female

    PWS.ASA.dat[[12]] <- as.matrix(PWS.ASA.dat[[12]][1:(nrow(PWS.ASA.dat[[12]])-years.to.remove),])    
    # MDM
    PWS.ASA.dat[[13]] <- as.matrix(PWS.ASA.dat[[13]][1:(nrow(PWS.ASA.dat[[13]])-years.to.remove),])    
    # egg deposition
    PWS.ASA.dat[[14]] <- as.matrix(PWS.ASA.dat[[14]][1:(nrow(PWS.ASA.dat[[14]])-years.to.remove),])    
    # egg deposition sd
    # 15 is first year of ADFG hydro survey. Needs no change.
    PWS.ASA.dat[[16]] <- as.matrix(PWS.ASA.dat[[16]][1:(nrow(PWS.ASA.dat[[16]])-years.to.remove),])    
    # ADFG hydroacoustic survey biomass
    # 17 is first year of PWSSC hydro survey. Needs no change.
    PWS.ASA.dat[[18]] <- as.matrix(PWS.ASA.dat[[18]][1:(nrow(PWS.ASA.dat[[18]])-years.to.remove),])    
    # PWSSC hydroacoustic survey biomass
    PWS.ASA.dat[[19]] <- as.matrix(PWS.ASA.dat[[19]][1:(nrow(PWS.ASA.dat[[19]])-years.to.remove),])    
    # PWSSC hydroacoustic survey biomass sd
    PWS.ASA.dat[[20]] <- PWS.ASA.dat[[20]][1:(nrow(PWS.ASA.dat[[20]])-years.to.remove),]               
    # seine age composition
    PWS.ASA.dat[[21]] <- PWS.ASA.dat[[21]][1:(nrow(PWS.ASA.dat[[21]])-years.to.remove),]               
    # spawner age composition
    PWS.ASA.dat[[22]] <- as.matrix(PWS.ASA.dat[[22]][1:(nrow(PWS.ASA.dat[[22]])-years.to.remove),])    
    # juvenile aerial survey

    write.data(PWS.ASA.dat, here::here(dir_retro, "PWS_ASA.dat"))

    # Write ESS data file with actual sample sizes (fixed).
    PWS.ASA.ESS.ctl <- dat.files$PWS_ASA_ESS.ctl 
    PWS.ASA.ESS.ctl[[2]] <- as.matrix(PWS.ASA.ESS.ctl[[2]][1:(nrow(PWS.ASA.ESS.ctl[[2]])-years.to.remove),])
    PWS.ASA.ESS.ctl[[3]] <- as.matrix(PWS.ASA.ESS.ctl[[3]][1:(nrow(PWS.ASA.ESS.ctl[[3]])-years.to.remove),])
    PWS.ASA.ESS.ctl[[4]] <- as.matrix(PWS.ASA.ESS.ctl[[4]][1:(nrow(PWS.ASA.ESS.ctl[[4]])-years.to.remove),])
    PWS.ASA.ESS.ctl[[5]] <- as.matrix(PWS.ASA.ESS.ctl[[5]][1:(nrow(PWS.ASA.ESS.ctl[[5]])-years.to.remove),])

    write.data(PWS.ASA.ESS.ctl, here::here(dir_retro, "PWS_ASA(ESS).ctl"))


    # Write covariate file (not currently used)
    PWS.ASA.covariate.ctl <- dat.files$PWS_ASA_covariate.ctl
    # Skip indices 1-4 for now, as they do not change with year
    PWS.ASA.covariate.ctl[[5]] <- as.matrix(PWS.ASA.covariate.ctl[[5]][1:(nrow(PWS.ASA.covariate.ctl[[5]])-years.to.remove),]) 
    PWS.ASA.covariate.ctl[[6]] <- as.matrix(PWS.ASA.covariate.ctl[[6]][1:(nrow(PWS.ASA.covariate.ctl[[6]])-years.to.remove),])
    # Skip indices 7-11 for now, as they do not change with year
    PWS.ASA.covariate.ctl[[12]] <- PWS.ASA.covariate.ctl[[12]][1:(nrow(PWS.ASA.covariate.ctl[[12]])-years.to.remove),]
    PWS.ASA.covariate.ctl[[14]] <- as.matrix(PWS.ASA.covariate.ctl[[14]][1:(nrow(PWS.ASA.covariate.ctl[[14]])-years.to.remove),])

    write.data(PWS.ASA.covariate.ctl, here::here(dir_retro, "PWS_ASA(covariate).ctl"))

    # Write age comp sample sizes file with actual 
    # sample sizes (currently fixed).
    age.comp.samp.sizes.txt <- dat.files$agecomp_samp_sizes.txt
    age.comp.samp.sizes.txt[[1]] <- as.matrix(age.comp.samp.sizes.txt[[1]][1:(nrow(age.comp.samp.sizes.txt[[1]])-years.to.remove),])
    age.comp.samp.sizes.txt[[2]] <- as.matrix(age.comp.samp.sizes.txt[[2]][1:(nrow(age.comp.samp.sizes.txt[[2]])-years.to.remove),])
    age.comp.samp.sizes.txt[[3]] <- as.matrix(age.comp.samp.sizes.txt[[3]][1:(nrow(age.comp.samp.sizes.txt[[3]])-years.to.remove),])
    age.comp.samp.sizes.txt[[4]] <- as.matrix(age.comp.samp.sizes.txt[[4]][1:(nrow(age.comp.samp.sizes.txt[[4]])-years.to.remove),])

    write.table(
        rbind(" ", " ", " ", " ", age.comp.samp.sizes.txt[[1]], " ", " "),
        file = here::here(dir_retro, "agecomp_samp_sizes.txt"), 
        append = F, sep = " ",
        row.names=FALSE, col.names=FALSE,
        quote=F
    )
    for(i in 2:length(age.comp.samp.sizes.txt)){
        write.table(
            rbind(age.comp.samp.sizes.txt[[i]], " ", " "),
            file = here::here(dir_retro, "agecomp_samp_sizes.txt"), 
            append = T, sep = " ",
            row.names=FALSE, col.names=FALSE,
            quote=F
        )
    }

    PWS.ASA.disease.dat <- dat.files$PWS_ASA_disease.dat
    PWS.ASA.disease.dat[[1]] <- PWS.ASA.disease.dat[[1]][1:(nrow(PWS.ASA.disease.dat[[1]])-years.to.remove),]  # Observed VHSV antibodies
    # Skip indices 2-4, as they do not change with year.
    PWS.ASA.disease.dat[[5]] <- PWS.ASA.disease.dat[[5]][1:(nrow(PWS.ASA.disease.dat[[5]])-years.to.remove),]   # Observed ichthyophonus antibodies
    # Skip indices 6-8, as they do not change with year.
    write.data(PWS.ASA.disease.dat, here::here(dir_retro, "PWS_ASA_disease.dat"))

    options(scipen=0)

}