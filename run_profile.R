#
# NOTE: This code only works for the seroprevalnce version
# of the BASA model. If using the non-sero model, you will
# need to modify lines 54-60. The processing code that 
# follows is not fully tested.
#

source(paste0(here::here("functions/"), "/run_basa.r"))

profile.dir <- file.path(here::here(), "profiles")
if(!dir.exists(profile.dir)){
    dir.create(profile.dir)
}

# this is the parameter names you woant to modify
# and needs to EXACTLY match the name a parameter
# in the PWS_ASA.par file. For now, it only works
# on single value parameters, not on vectors.
par.name <- "Z_0_8"

profile.values <- seq(0.1, 0.8, 0.025)

setwd(profile.dir)

profile.par.dir <- file.path(profile.dir, par.name)
if(!dir.exists(profile.par.dir)){
    dir.create(profile.par.dir)
}

overwrite <- TRUE

for(v in profile.values){
    dir.name <- paste0(par.name, "=", v)
    profile.value.dir <- file.path(profile.par.dir, dir.name)
    if(!dir.exists(profile.value.dir)){
        dir.create(profile.value.dir)
    }

    if(file.exists(file.path(profile.value.dir, "PWS_ASA.rep")) & !overwrite){
        print(paste("Skipped", par.name, "=", v))
        next;
    }

    file.list <- c("PWS_ASA.tpl", "PWS_ASA.cor", "PWS_ASA.PIN", "PWS_ASA.par", "agecomp_samp_sizes.txt", "PWS_ASA.dat", "PWS_ASA(covariate).ctl", "PWS_ASA(ESS).ctl", "PWS_ASA(ESS_estimate).ctl", "PWS_ASA(phases).ctl", "PWS_ASA(sim_settings).ctl", "PWS_ASA_disease.dat", "PWS_ASA(par).ctl")
    file.copy(
        from = file.path(here::here(), "model", file.list),
        to = profile.value.dir,
        overwrite = TRUE,
        copy.mode = TRUE
    )

    setwd(profile.value.dir)

    lines <- readLines("PWS_ASA.tpl")
    #idxs <- grepl("_PIN = ", lines)
    #PIN.lines <- lines[idxs]
    par.grepl <- grepl(paste0(par.name, "_PIN = "), lines)
    lines[par.grepl] <- gsub("(?<= = )[\\d\\.]+", v, lines[par.grepl], perl=TRUE)

    writeLines(lines, "PWS_ASA.tpl")

    ##################
    ## Just need to run in MLE mode to get likelihood values
    ##################
    system("admb -s PWS_ASA")
    system("./PWS_ASA -pinwrite")

    setwd(profile.dir)

}

# profile.dirs <- file.path(here::here(), "profiles", par.name, paste0(par.name, "=", profile.values))
# likelihoods <- file.path(profile.dirs, "mcmc_out", "llikcomponents.csv")

# base.model.liketable <- read_csv(file.path(here::here(), "model", "mcmc_out", "llikcomponents.csv"), col_names = FALSE) %>%
#     `colnames<-`(c("SeAC_Like", "SpAC_Like", "Egg_Like", "ADFGhyd_Like", "PWSSChyd_Like", "MDM_Like", "age0pen_Like", "mortdevspen_Like", "age0covar_prior_Like", "mortcovar_prior_Like", "Z_prior_Like", "ADFGhyd_prior_Like", "PWSSChyd_prior_Like", "m_prior_Like", "juv_Like", "sero_Like", "full_Like")) %>%
#     mutate(
#         full = full_Like,
#         agecomp = -SeAC_Like+-SpAC_Like,
#         survey = Egg_Like+ADFGhyd_Like+PWSSChyd_Like+MDM_Like+juv_Like,
#         seroprev = sero_Like,
#         pen = age0pen_Like+mortdevspen_Like,
#         prior = age0covar_prior_Like+mortcovar_prior_Like+Z_prior_Like+ADFGhyd_prior_Like+PWSSChyd_prior_Like+m_prior_Like,
#     ) %>%
#     select(-c("SeAC_Like", "SpAC_Like", "Egg_Like", "ADFGhyd_Like", "PWSSChyd_Like", "MDM_Like", "age0pen_Like", "mortdevspen_Like", "age0covar_prior_Like", "mortcovar_prior_Like", "Z_prior_Like", "ADFGhyd_prior_Like", "PWSSChyd_prior_Like", "m_prior_Like", "juv_Like", "sero_Like", "full_Like")) %>%
#     pivot_longer(everything(), names_to="component", values_to="LLike")


# make.likelihood.tibble <- function(fname, par.value){
#     return(
#         read_csv(fname, col_names = FALSE) %>%
#             `colnames<-`(c("SeAC_Like", "SpAC_Like", "Egg_Like", "ADFGhyd_Like", "PWSSChyd_Like", "MDM_Like", "age0pen_Like", "mortdevspen_Like", "age0covar_prior_Like", "mortcovar_prior_Like", "Z_prior_Like", "ADFGhyd_prior_Like", "PWSSChyd_prior_Like", "m_prior_Like", "juv_Like", "sero_Like", "full_Like")) %>%
#             mutate(
#                  full = full_Like,
#                  agecomp = -SeAC_Like+-SpAC_Like,
#                  survey = Egg_Like+ADFGhyd_Like+PWSSChyd_Like+MDM_Like+juv_Like,
#                  seroprev = sero_Like,
#                  pen = age0pen_Like+mortdevspen_Like,
#                  prior = age0covar_prior_Like+mortcovar_prior_Like+Z_prior_Like+ADFGhyd_prior_Like+PWSSChyd_prior_Like+m_prior_Like
#             ) %>%
#             select(-c("SeAC_Like", "SpAC_Like", "Egg_Like", "ADFGhyd_Like", "PWSSChyd_Like", "MDM_Like", "age0pen_Like", "mortdevspen_Like", "age0covar_prior_Like", "mortcovar_prior_Like", "Z_prior_Like", "ADFGhyd_prior_Like", "PWSSChyd_prior_Like", "m_prior_Like", "juv_Like", "sero_Like", "full_Like")) %>%
#             pivot_longer(everything(), names_to="component", values_to="LLike") %>%
#             mutate(
#                 base.LLike = base.model.liketable %>% pull(LLike) 
#             ) %>%
#             mutate(
#                 LLike.diff = LLike - base.LLike,
#                 par.val = par.value
#             ) %>%
#             select(-c(base.LLike))
#     )
# }

# likelihood.tibble <- tibble()
# for(i in 1:length(likelihoods)){
#     fname <- likelihoods[i]
#     val <- profile.values[i]

#     likelihood.tibble <- likelihood.tibble %>% bind_rows(make.likelihood.tibble(fname, val))

# }

# lt <- likelihood.tibble %>%
#         mutate(
#             component = factor(component, levels=c("full", "agecomp", "survey", "seroprev", "prior", "pen"), labels=c("Full", "Age Composition", "Surveys", "Seroprevalence", "Priors", "Penalties"))
#         ) %>%
#         group_by(component, par.val) %>%
#         median_qi(LLike.diff, .width=c(0.50, 0.95)) %>%
#         print(n=100)

# ggplot(lt)+
#     geom_line(aes(x=par.val, y=LLike.diff, color=component), size=1)+
#     geom_point(aes(x=par.val, y=LLike.diff, color=component), size=3)+
#     geom_hline(aes(yintercept=-2), linetype="dashed")+
#     geom_hline(aes(yintercept=2), linetype="dashed")+
#     geom_hline(aes(yintercept=0))+
#     #scale_color_manual(values=c("black", "red", "blue", "#009500", "grey60", "grey30")) + 
#     #scale_y_continuous(limits=c(-5, 5))+
#     coord_cartesian(ylim=c(-10, 10))+
#     theme_classic()

#load(file.path(here::here(), "profiles", "Z_0_8", "Z_0_8=0.15", "mcmc_out", "NUTS_fit.RDS"))

par.name <- "Z_0_8"
profile.values <- seq(0.1, 0.775, 0.025)
profile.dirs <- file.path(here::here(), "profiles", par.name, paste0(par.name, "=", profile.values))

rep <- readLines(file.path(here::here(), "model", "rep_out", "PWS_ASA.rep"), n=40)
likes <- rep[grepl("[[:<:]][\\d\\.]+", rep, perl=TRUE)]
likes <- as.numeric(likes[2:length(likes)])
names(likes) <- c("Full", "SeAC", "SpAC", "Egg", "ADFGhyd", "PWSSChyd", "mdm", "age0_pen", "mort_pen", "age0_prior", "mort_prior", "Z_prior", "ADFGhyd_prior", "PWSSChyd_prior", "m_prior", "juv", "sero")
base.likes <- as_tibble(data.frame(as.list(likes))) %>%
    mutate(
        agecomp = SeAC+SpAC,
        survey = Egg+ADFGhyd+PWSSChyd+mdm+juv,
        seroprev = sero,
        pen = age0_pen+mort_pen,
        prior = age0_prior+mort_prior+Z_prior+ADFGhyd_prior+PWSSChyd_prior+m_prior
    ) %>%
    select(Full, agecomp, survey, seroprev, pen, prior) %>%
    pivot_longer(everything(), names_to="component", values_to="Like")

rep.files <- file.path(profile.dirs, "rep_out", "PWS_ASA.rep")

likelihood.df <- data.frame()
biomass.data <- NA
for(f in rep.files){
    rep <- readLines(f)
    rep.1 <- rep[1:40]
    likes <- rep.1[grepl("[[:<:]][\\d\\.]+", rep.1, perl=TRUE)]
    likes <- as.numeric(likes[2:length(likes)])
    names(likes) <- c("Full", "SeAC", "SpAC", "Egg", "ADFGhyd", "PWSSChyd", "mdm", "age0_pen", "mort_pen", "age0_prior", "mort_prior", "Z_prior", "ADFGhyd_prior", "PWSSChyd_prior", "m_prior", "juv", "sero")

    likelihood.df <- bind_rows(likelihood.df, data.frame(as.list(likes)))

    biomass <- as.numeric(strsplit(rep[149], " ")[[1]])
    final.biomass <- biomass[length(biomass)]
    biomass.data <- c(biomass.data, final.biomass)
}

piner <- as_tibble(likelihood.df) %>%
    mutate(
        agecomp = SeAC+SpAC,
        survey = Egg+ADFGhyd+PWSSChyd+mdm+juv,
        seroprev = sero,
        pen = age0_pen+mort_pen,
        prior = age0_prior+mort_prior+Z_prior+ADFGhyd_prior+PWSSChyd_prior+m_prior
    ) %>%
    select(Full, agecomp, survey, seroprev, pen, prior) %>%
    pivot_longer(everything(), names_to="component", values_to="Like") %>%
    mutate(
        par.val = rep(profile.values, each=6)
    ) %>%
    left_join(base.likes, by="component") %>%
    mutate(LIKE = Like.x - Like.y) %>%
    select(-c(Like.x, Like.y)) %>%
    mutate(
        component=factor(component, levels=c("Full", "agecomp", "survey", "seroprev", "prior", "pen"), labels=c("Full", "Age Compositions", "Surveys", "Seroprevalence", "Priors", "Penalties"))
    )

ggplot(piner)+
    geom_line(aes(x=par.val, y=LIKE, color=component), size=1)+
    geom_point(aes(x=par.val, y=LIKE, color=component), size=3)+
    geom_hline(aes(yintercept=-2), linetype="dashed")+
    geom_hline(aes(yintercept=2), linetype="dashed")+
    geom_hline(aes(yintercept=0))+
    #scale_color_manual(values=c("black", "red", "blue", "#009500", "grey60", "grey30")) + 
    #scale_y_continuous(limits=c(-5, 5))+
    coord_cartesian(ylim=c(-100, 100))+
    labs(x=par.name, y="Change in Likelihood", color="Likelihood \nComponent", title=paste("Likelihood Profile over", par.name))+
    theme_classic()+
    theme(
        axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        plot.title = element_text(size=18)
    )

as_tibble(biomass.data) %>% na.omit() %>%
    mutate(par.val = profile.values) %>%

    ggplot()+
        geom_line(aes(x=par.val, y=value))+
        geom_hline(aes(yintercept=20000), linetype="dashed")+
        geom_hline(aes(yintercept=40000), linetype="dashed")+
        scale_y_continuous(labels = scales::comma)+
        coord_cartesian(ylim=c(0, 40000))+
        theme_classic()
