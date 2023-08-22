source(file=file.path(here::here(), "functions", "fun_read_dat.R"))
source(file=file.path(here::here(), "functions", "simulation", "operating_model", "fun_fish.R"))
source(file=file.path(here::here(), "functions", "simulation", "operating_model", "fun_operm.R"))
files.sources = list.files(here::here("functions/simulation/operating_model/control_rules"), full.names=TRUE)
sapply(files.sources, source)

nage=10

listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}

initialize.popdyn.variables <- function(nyr.sim){

  rownames <- 1:nyr.sim
  colnames <- 0:(nage-1)

  survival.summer         <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  survival.winter         <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  maturity                <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  prefish.spawn.biomass   <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  seine.catch             <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  gillnet.catch           <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  pound.catch             <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  foodbait.catch          <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  n.spawners              <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  spawn.biomass.age.comp  <- matrix(0, nyr.sim, nage,   dimnames=list(rownames, colnames))
  spawn.biomass           <- rep(0, nyr.sim)
  seine.biomass           <- rep(0, nyr.sim)
  true.nya                <- matrix(0, nyr.sim+1, nage, dimnames=list(c(rownames, nyr.sim+1), colnames))
  annual.age0.devs        <- rep(0, nyr.sim)

  pop_dyn <- listN(survival.summer, survival.winter, maturity, prefish.spawn.biomass, 
                seine.catch, gillnet.catch, pound.catch, foodbait.catch,
                n.spawners, spawn.biomass.age.comp, true.nya,
                annual.age0.devs)

  return(pop_dyn)
  
}

set.initial.conditions <- function(dir, pop_dyn, sim.seed, start.year=2023){ 
  # Takes a random sample of the NYA distribution from year 0
  # (either then actual stock assessment or the hindcast model).
  set.seed(sim.seed)
  nyr <- start.year-1980+1
  init.nya <- read_csv(paste0(dir, "/mcmc_out/Num_at_age.csv"), col_names=FALSE, show_col_types = FALSE) %>%
                select_at(((10*nyr)-9):(10*nyr)) %>%
                slice_sample(n=1)

  pop_dyn$true.nya[1, ] <- as.numeric(init.nya) # Should these be rounded to integers?
  # pop_dyn$prefish.spawn.biomass[1, ] <- mat*as.numeric(init.nya)*waa
  
  return(pop_dyn)
}



run.simulation <- function(hcr.options, nyr.sim, sim.seed=NA, write=NA, 
                                  start.year=1, stop.year=NA, init.start.year=2021,
                                  assessment=TRUE, hindcast=TRUE, cr.name="test", max.regime.length=15,
                                  custom.recruitment.devs = NA){
    # print(write)
    if(is.na(write)){
      write <- paste0(here::here("results"), "/test")
    }

    # Copy over current stock assessment to use for starting values
    basa.root.dir <- "~/Desktop/Projects/basa/model/"
    model.0.dir <- basa.root.dir
    write <- model.0.dir

    # Read in data files JUST ONCE, store then write within code
    dat.files <- read.data.files(model.0.dir)

    # Read these out of the mcmc_out/*.csv files
    waa.data <- dat.files$PWS_ASA.dat[[4]]
    true.waa <- apply(waa.data[31:41, ], 2, median)

    fec.data <- dat.files$PWS_ASA.dat[[5]]
    true.fec <- apply(fec.data, 2, median, na.rm=TRUE)
    
    # Assuming knife-edged selectivity at age-3
    fish.selectivity <- matrix(1, nrow=4, ncol=10)
    fish.selectivity[, 1:3] <- 0 # Selectivity 0 for fish age 0-2

    log.mean.age0 <- read_csv(file.path(model.0.dir, "mcmc_out", "Num_at_age.csv"), col_names = FALSE, show_col_types = FALSE) %>%
      select(seq(1, 430, 10)) %>%
      rename_with(~ as.character(1980:2022), everything()) %>%
      summarise(
        across(
          as.character(1980:2022), 
          \(x) mean(x, na.rm=TRUE)
        )
      ) %>%
      as.matrix %>%
      as.vector %>%
      mean %>%
      log

    # Read in other parameters
    params <- read.par.file("~/Desktop/Projects/basa/model/PWS_ASA.par")
    par.samples <- read_csv(paste0(model.0.dir, "/mcmc_out/iterations.csv"), show_col_types=FALSE) %>% 
                na.omit() %>%
                select(!matches("[0-9]]$")) %>%
                slice_sample(n=1) %>% as.list

    params[names(par.samples)] <- par.samples
    params$log_MeanAge0     <- log.mean.age0
    params$female.spawners  <- tail(dat.files$PWS_ASA.dat$perc.female, 1)
    params$pk               <- dat.files$PWS_ASA.dat$pk
    params$waa              <- true.waa
    params$fec              <- true.fec
    params$selectivity      <- fish.selectivity
    params$catch.sd         <- 0.1

    maturity <- calc.maturity(params$mat_par_1, params$mat_par_2)

    # Initialize projection estimates (e.g. from forecast model)
    pop_dyn <- initialize.popdyn.variables(nyr.sim)
    pop_dyn <- set.initial.conditions(write, pop_dyn, sim.seed, init.start.year-1)

    # Generate new age-0 recruitment deviates.
    # Devs are pulled from a regime-based recruitment function.
    recruitment <- custom.recruitment.devs
    print(recruitment)
    pop_dyn$annual.age0.devs <- recruitment$devs
    params$sigma_age0devs <- recruitment$sigmas
    #pop_dyn$annual.age0.devs <- rnorm(nyr.sim, mean=-0.35, sd=0.9) # change this so that the fishery doesn't immediately recover
    print(pop_dyn$annual.age0.devs)

    # Start loop
    control.rule <- rep(0, nyr.sim)

    if(is.na(stop.year)){
      stop.year <- nyr.sim
    }

    for(y in start.year:stop.year){  

        true.ssb <- sum(maturity*pop_dyn$true.nya[y, ]*params$waa)
        true.nya <- pop_dyn$true.nya[y, ]

        # If we want to run the full assessment, then use the assessment estimates in the HCR
        # calculation. Otherwise, use the deterministic biomass. This is useful for debugging
        # the operating model quickly (no need to run the whole assessment).
        hcr.name <- hcr.options$type

        hcr.params <- list(
              curr.biomass = true.ssb,
              age.structure = true.nya,
              rel.biomass = true.ssb/sum(pop_dyn$prefish.spawn.biomass[y-3, ]),
              options = hcr.options
            )

        tar_hr <- do.call(hcr.name, c(hcr.params))

        control.rule[y] <- tar_hr

        # Execute fishery following management recommendation
        catch.at.age <- fun_fish(tar_hr, true.ssb, true.nya, params$waa, params$selectivity, params$catch.sd)

        # Run operating model
        pop_dyn <- fun_operm(y, pop_dyn$true.nya[y, ], catch.at.age, params, pop_dyn, sim.seed, project=TRUE)

        # Read in the new data files and update the female.spawners param
        #dat.files <- read.data.files(model.dir)
        params$female.spawners  <- tail(dat.files$PWS_ASA.dat$perc.female, 1)

        print(paste0(cr.name, " (sim ", sim.seed, "): ", y, "/", stop.year))

    }
    return(list(pop.dyn=pop_dyn, harvest.rate=control.rule))
}