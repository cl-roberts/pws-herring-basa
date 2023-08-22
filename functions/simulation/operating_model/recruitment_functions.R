custom.recruitment.1 <- function(nyr.sim, sim.seed, max.regime.length=15){
    
  set.seed(sim.seed)
  basa.rec.devs <- read.par.file(file.path(here::here(), "model", "PWS_ASA.par"))$annual_age0devs
  return(list(
    devs=sample(basa.rec.devs, nyr.sim, replace=TRUE),
    sigmas=rep(0, nyr.sim)
  ))

}

custom.recruitment.2 <- function(nyr.sim, sim.seed, max.regime.length=15){

  set.seed(sim.seed)
  basa.rec.devs <- read.par.file(file.path(here::here(), "model", "PWS_ASA.par"))$annual_age0devs
  
  high <- list(mu=mean(basa.rec.devs[1:12]), sd=sd(basa.rec.devs[1:12]))
  low  <- list(mu=mean(basa.rec.devs[13:35]), sd=sd(basa.rec.devs[13:35]))

  regime <- sample(0:1, 1)

  devs <- rep(NA, nyr.sim)
  sigmas <- rep(NA, nyr.sim)

  # 1 --> High regime
  for(y in 1:(nyr.sim)){
    if(y %% max.regime.length == 0) regime <- !regime
    dev <- ifelse(regime == 1, rnorm(1, high$mu, high$sd), rnorm(1, low$mu, low$sd))
    sig <- ifelse(regime == 1, high$sd, low$sd)
    devs[y] <- dev
    sigmas[y] <- sig
  }

  return(list(
    devs=devs,
    sigmas=sigmas
  ))

}

custom.recruitment.3 <- function(nyr.sim, sim.seed, max.regime.length=15){
  set.seed(sim.seed)
  devs <- rep(NA, nyr.sim)
  sigmas <- rep(NA, nyr.sim)
  #devs[1:max.regime.length] <- rnorm(max.regime.length, 0.0645, 1.35)
  #sigmas[1:max.regime.length] <- rep(1.20, max.regime.length)
  high.regime <- sample(0:1, 1)
  for(y in 1:(nyr.sim)){
    if(y %% max.regime.length == 0) high.regime <- !high.regime
    dev <- ifelse(high.regime == 1, rnorm(1, -1.270, 1.05), rnorm(1, 0.0645, 1.35))
    sig <- ifelse(high.regime == 1, 1.05, 1.35)
    devs[y] <- dev
    sigmas[y] <- sig
  }
  return(list(devs=devs, sigmas=sigmas))
}
