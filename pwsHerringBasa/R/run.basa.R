#' Run BASA model
#'
#' Compiles ADMB code, reads in data, calculates effective sample sizes, and
#' executes BASA model using NUTS sampling algorithm. This function *should*
#' work for linux and mac, but has only been tested on windows.
#'
#' @param dir_model Relative path to model
#' @param n.samples Number of samples to draw from MCMC iterations
#' @param n.warmup Number of MCMC burn-in iterations
#' @param n.time Number of minutes after which the model will quit running
#' @param n.chains Number of MCMC chains to run
#'
#' @returns Metadata from model fit in ADMB. Model outputs are saved to `model/mcmc_out/`
#' and `model/rep_out` subdirectories.
#'

run.basa <- function(dir_model, n.samples=2000, n.warmup=700, n.time=5, n.chains=4){
 
  template.files <- here::here(dir_model)
  print(template.files)
  setwd(template.files)
  system("admb -s PWS_ASA")
  #shell("admb -s PWS_ASA")

  # identify operating system for shell commands
  OS <- switch(Sys.info()[['sysname']],
               Windows = "PC",
               Darwin = "Mac",
               Linux = "Linux")

  model.data <- read.admb.files()

  nyr.fit <- model.data$nyr.fit # 43 for data up through 2022
  nyr.tot <- model.data$nyr.tot # 43 for data up through 2022
  nage <- 10                    # number of age classes (age 0 to 9+) (only care about 3-9+)

  # Read in measured age comps
  seine.age.comp <- model.data$seine.age.comp
  spawn.age.comp <- model.data$spawn.age.comp
  vhsv.age.comp <- model.data$vhsv.age.comp
  ich.age.comp <- model.data$ich.age.comp

  # Read in the actual sample sizes
  seine.samp.size <- read.data.files(template.files)$"agecomp_samp_sizes.txt"$seine_sample_size
  spawn.samp.size <- read.data.files(template.files)$"agecomp_samp_sizes.txt"$spawn_sample_size
  vhsv.samp.size  <- read.data.files(template.files)$"agecomp_samp_sizes.txt"$vhsv_sample_size
  ich.samp.size   <- read.data.files(template.files)$"agecomp_samp_sizes.txt"$ich_sample_size

  # The following commands skips over lines of the file to read in tables
  # so change the number skipped if file is modified
  # seine.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4,                 nrows=nyr.tot)
  # spawn.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+1*(nyr.tot+1),   nrows=nyr.tot)
  # vhsv.samp.size  <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+2*(nyr.tot+1),   nrows=nyr.tot)
  # ich.samp.size   <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+3*(nyr.tot+1)+1, nrows=nyr.tot)


  # Create empty matrices to fill estimated ESS and age comps
  seine.ess.its <- matrix(0, nyr.tot, 1)  # Matrix to hold all iterations of the routine
  seine.ess.its <- seine.samp.size        # fill in the first column with the recorded sample size

  spawn.ess.its <- matrix(0, nyr.tot, 1)
  spawn.ess.its <- spawn.samp.size

  vhsv.ess.its  <- matrix(0, nyr.tot, 1)
  vhsv.ess.its  <- vhsv.samp.size
  ich.ess.its   <- matrix(0, nyr.tot, 1)
  ich.ess.its   <- ich.samp.size

  # Change phases of the ESS in phases file to use the PWS_ASA(ESS_estimate)
  phases <- readLines("PWS_ASA(ESS).ctl", -1)
  phases[5] <- 1
  writeLines(phases, "PWS_ASA(ESS).ctl")

  # LOOP THROUGH AND ITERATIVELY CALCULATE ESS
  convergence <- 0

  seine.ess <- seine.samp.size
  spawn.ess <- spawn.samp.size
  vhsv.ess  <- vhsv.samp.size
  ich.ess   <- ich.samp.size
  its <- 1

  age.comps <- list(seine=seine.age.comp, spawn=spawn.age.comp, vshv=vhsv.age.comp, ich=ich.age.comp)
  start.ess <- list(seine=seine.ess, spawn=spawn.ess, vhvs=vhsv.ess, ich=ich.ess)
  samp.size <- list(seine=seine.samp.size, spawn=spawn.samp.size, vhsv=vhsv.samp.size, ich=ich.samp.size)

  # calc.ess <- calculate.ess(age.comps, start.ess, samp.size, nyr.fit)

  for(i in 1:2){
    # Create "PWS_ASA(ESS_estimate).ctl" with sample sizes (the original sample sizes on the first iteration)
    write.table(
      rbind("# PWS age comp effective sample sizes",
            "# Seine ESS",      seine.ess,  " ",
            "# Spawn ESS",      spawn.ess,  " ",
            "# VHS sero ESS",   vhsv.ess,   " ",
            "# Ich prev ESS",   ich.ess
      ),
      file = "PWS_ASA(ESS_estimate).ctl",
      append = F,
      sep = " ",
      row.names=FALSE,
      col.names=FALSE,
      quote=F
    )

    # Compile and Run PWS_ASA
    if((OS == "Mac") | (OS == "Linux")){
      # system("./PWS_ASA -pinwrite -nohess")
      system("./PWS_ASA -pinwrite")
    }else if(OS == "PC"){
      shell("PWS_ASA  -pinwrite -nohess")
    }


    # Read in the estimated seine and spawner age comps
    seine.age.comp.est <- read.table("rep_out/SeAC_pd.rep",     header = FALSE)
    spawn.age.comp.est <- read.table("rep_out/SpAC_pd.rep",     header = FALSE)
    vhsv.age.comp.est  <- read.table("rep_out/vhscomp_pd.rep",  header = FALSE)
    ich.age.comp.est   <- read.table("rep_out/ichcomp_pd.rep",  header = FALSE)

    # Calculate the ESS
    seine.ess <- rowSums(seine.age.comp.est*(1-seine.age.comp.est))/rowSums((seine.age.comp[1:nyr.fit, ]-seine.age.comp.est)^2)
    spawn.ess <- rowSums(spawn.age.comp.est*(1-spawn.age.comp.est))/rowSums((spawn.age.comp[1:nyr.fit, ]-spawn.age.comp.est)^2)
    vhsv.ess  <- rowSums(vhsv.age.comp.est*(1-vhsv.age.comp.est))/rowSums((vhsv.age.comp[1:nyr.fit, ]-vhsv.age.comp.est)^2)
    ich.ess   <- rowSums(ich.age.comp.est*(1-ich.age.comp.est))/rowSums((ich.age.comp[1:nyr.fit, ]-ich.age.comp.est)^2)

    # Remove the missing years of age comps
    seine.ess.rem <- seine.ess[!(seine.age.comp[1:nyr.fit, 1]==-9)]
    spawn.ess.rem <- spawn.ess[!(spawn.age.comp[1:nyr.fit, 1]==-9)]
    vhsv.ess.rem <- vhsv.ess[!(vhsv.age.comp[1:nyr.fit, 1]==-9)]
    ich.ess.rem <- ich.ess[!(ich.age.comp[1:nyr.fit, 1]==-9)]

    # Calculate the ratio of ESS to original sample sizes
    seine.ratio <- seine.ess.rem/seine.samp.size[1:nyr.fit, 1][seine.age.comp[1:nyr.fit, 1]!=-9]
    spawn.ratio <- spawn.ess.rem/spawn.samp.size[1:nyr.fit, 1][spawn.age.comp[1:nyr.fit, 1]!=-9]
    vhsv.ratio  <- vhsv.ess.rem/vhsv.samp.size[1:nyr.fit, 1][vhsv.age.comp[1:nyr.fit, 1]!=-9]
    ich.ratio   <- ich.ess.rem/ich.samp.size[1:nyr.fit, 1][ich.age.comp[1:nyr.fit, 1]!=-9]

    # Calculate the harmonic means (see Muradian et al. 2017 and Stewart & Hamel 2014)
    seine.hm <- 1/mean(1/seine.ratio)
    spawn.hm <- 1/mean(1/spawn.ratio)
    vhsv.hm <- 1/mean(1/vhsv.ratio)
    ich.hm <- 1/mean(1/ich.ratio)

    # Compare this harmonic mean to the previous using a convergence criteria (WHAT AM I CONVERGING!!!!)
    if(its==1) {
      convergence <- 0
      seine.hmS <- seine.hm
      spawn.hmS <- spawn.hm
      vhsv.hmS  <- vhsv.hm
      ich.hmS   <- ich.hm
    }else{
      seine.test <- abs(seine.hm - seine.hmS[its-1])/seine.hmS[its-1]*100
      spawn.test <- abs(spawn.hm - spawn.hmS[its-1])/spawn.hmS[its-1]*100
      vhsv.test  <- abs(vhsv.hm - vhsv.hmS[its-1])/vhsv.hmS[its-1]*100
      ich.test   <- abs(ich.hm - ich.hmS[its-1])/ich.hmS[its-1]*100

      convergence <- (seine.test<0.1 & spawn.test<0.1) # This criteria was arbitrarily chosen (0.1% change)
      seine.hmS <- rbind(seine.hmS, seine.hm)
      spawn.hmS <- rbind(spawn.hmS, spawn.hm)
      vhsv.hmS  <- rbind(vhsv.hmS, vhsv.hm)
      ich.hmS   <- rbind(ich.hmS, ich.hm)
    }

    # Now multiply the harmonic mean by the sample size to get the new ESS
    seine.ess <- round(seine.hm*seine.samp.size, 0)
    spawn.ess <- round(spawn.hm*spawn.samp.size, 0)
    vhsv.ess  <- round(vhsv.hm*vhsv.samp.size,   0)
    ich.ess   <- round(ich.hm*ich.samp.size,     0)

    # Use the average ESS for all years (each years obs weighted equally)
    # seine.ess[seine.ess>0] <- round(mean(seine.ess[seine.ess>0]), digits=0)
    # spawn.ess[spawn.ess>0] <- round(mean(spawn.ess[spawn.ess>0]), digits=0)

    # Denote the missing values
    seine.ess[(seine.age.comp[, 1] == -9), 1] <- -9
    spawn.ess[(spawn.age.comp[, 1] == -9), 1] <- -9
    vhsv.ess[(vhsv.age.comp[, 1]   == -9), 1] <- -9
    ich.ess[(ich.age.comp[, 1]     == -9), 1] <- -9

    vhsv.ess <- vhsv.samp.size
    ich.ess <- ich.samp.size

    # Fill in this iteration"s ESS
    seine.ess.its <- cbind(seine.ess.its, round(seine.ess, 0))
    spawn.ess.its <- cbind(spawn.ess.its, round(spawn.ess, 0))
    vhsv.ess.its <- cbind(vhsv.ess.its, round(vhsv.ess, 0))
    ich.ess.its <- cbind(ich.ess.its, round(ich.ess, 0))

    # Cease iterations if convergence hasn"t happened after so many...
    if(its==10) {
      break
    }
    its <- its+1
  }

  # Turn of the phases for the ESS calculation so it no longer recalculates ESS in future model runs
  # Now write the converged ESS to a ctl file to be used for model runs
  write.table(
    rbind(
      "# PWS age comp effective sample sizes", paste0("# (", date(), ")"), " ",
      "# Determines which ctl file for the age comps and ESS to use (1 uses ESS control to be iteratively estimated)", -1, " ",
      "# Seine ESS", seine.ess, " ",
      "# Spawn ESS", spawn.ess, " ",
      "# Sero ESS", vhsv.ess, " ",
      "# Ich prev ESS", ich.ess
    ),
    file = "PWS_ASA(ESS).ctl",
    append = F,
    sep = " ",
    row.names=FALSE,
    col.names=FALSE,
    quote=F
  )

  ######################################################
  # Create reps x starting par vectors, and run NUTS
  setwd(template.files)
  reps <- n.chains
  set.seed(8558)
  seeds <- sample(1:1e4, size=reps)
  #system("admb -s PWS_ASA")
  system("./PWS_ASA -pinwrite -hbf 1")

  inits <- init.admb.params(reps)

  # ADMB command for running nuts
  # PWS_ASA -nox -noest -nohess -maxfn 0 -nuts -mcmc 2000 -warmup 500 -chain 1 -mcseed 8682524 -max_treedepth 12 -adapt_delta 0.8 -adapt_mass -mcpin init.pin

  # Pilot run to check

  if(dir.exists("mcmc_out")){
    system("rm mcmc_out/*.csv")
  }

  start.time <- Sys.time()
  fit.1 <- adnuts::sample_nuts(model='./PWS_ASA',path=template.files,
                       iter=n.samples,
                       warmup=n.warmup,
                       #warmup=100,
                       duration = n.time,
                       init=inits, seeds=seeds, chains=reps,cores=reps,
                       mceval=TRUE,
                       control=list(
                         adapt_delta=0.9,
                         #max_treedepth=16,
                         metric="mle"
                       )
  )
  end.time <- Sys.time()
  total.time <- end.time - start.time

  return(list(fit1=fit.1, time=total.time))
}
