
# changelog

## CLR - 6/30/2025

- Fixed mortality covariate indexing and mortality forecast bugs in ADMB model

## CLR - 7/1/2025

- Fixed bug in estimation model for juvenile schools, now better fits
- Estimate `log_MeanAge0`, previously a fixed parameter (as is in ADMB model)
- Estimate age0 deviate for last year in model instead of fixing to 0
- Made age0 deviates random effects
- scaled age0 deviates by sigma_age0devs in attempt to decorrelate deviates
 
## CLR - 7/2/2025

- Added disease prevalence data from 2020-2024

## CLR - 7/3/2025

- undo treating age0 deviates as random effects in TMB model
- re-fix last deviate to 0


## CLR - 7/5/2025

- fixed divergences related to posfun in TMB model
- TMB bug fixes

## CLR - 7/6/2025

- TMB bug fixes

## CLR - 7/7/2025

- TMB bug fixes

## CLR - 7/9/2025

- implemented ESS calculation in TMB model
- implement report writing function in TMB model for debugging

## CLR - 7/10/2025

- fixed pound spawn removals calculation in first year of ADMB model (for calculating post-fishery spawning NAA) 
- fixed ignoring 7+ catches in first year of ADMB model
- fixed plus group bug (first five years) in TMB model 


## CLR - 7/11/2025

- fixed incrementation in first five years of age comp likelihood calculations in both ADMB and TMB models


## CLR - 7/23/2025

- fixed likelihood penalties applied by posfun in ADMB model

## CLR - 8/9/2025

- move TMB model to main
- TMB model is now the main model

## CLR - 8/12/2025

- update documentation