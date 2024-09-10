################################################################################

# Script for performing retrospective analyses on BASA

# Creates a new subdirectory ('retrospectives/') and runs
# BASA on successively fewer years of data (up to a provided
# number of 'peels'). Final biomass estimates and recruitment
# deviates are aggregated and analysed for bias patterns relative
# to the most recent (current year) full BASA run.

# authors: Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and max number of years to peel 
#         (n.peels)

# outputs: BASA outputs for each retrospective peel. Each peeled run is saved to a 
#          subdirectory of 'retrospectives/'.

################################################################################

# attach packages ----

library(doParallel)
library(pwsHerringBasa)

# directory handling ----

dir_model <- here::here("model")
dir_retro <- here::here("retrospectives")

# local variables ----

n.peels <- 5
run.parallel <- FALSE
total.cores <- parallel::detectCores()
parallel.runs <- (total.cores-1) %/% 4      # BASA runs using 4 nodes

# run retrospective analysis ----

# Need to copy and then modify all of input files so they only reflect the data
# that was available in that year. Then actually run the assessment on all of
# the old datasets. Do this in parallel to speed things up.
if(run.parallel){

    # CLR: NOT TESTED ---- 
    cluster <- parallel::makeCluster(parallel.runs, type="FORK")
    registerDoParallel(cluster)
    foreach(i=1:n.peels) %dopar% {
        run.basa.retro(i, dir_model, dir_retro)
    }
    stopCluster(cluster)

} else {

    for(i in 1:n.peels){
        run.basa.retro(i, dir_model, dir_retro)
    }

}
