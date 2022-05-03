# Prince William Sound Herring Bayesian Age-Structured Stock Assessment Model (BASA)

The Bayesian Age-Structured Stock Assessment Model (BASA) for Pacific herring is the succesor to the historical Age-Structured Assessment model (ASA) used by the Alaska Department of Fish and Game (ADF&G) to assess the status of Pacific herring (_Clupea pallasii_) in Prince William Sound. Historically, the PWS herring stock posessed a total stock biomass in upwards of 100,000 metric tons, more than anough to sustain several active fisheries. The stock suddenly crashed in the spring of 1993, resulting in a complete closure of the fishery. To this day, the stock biomass remains below 20,000 metric tons, preventing the opening of the fishery per the harvest control rule put in place by ADF&G in 1994.

BASA has been under active development by graduate students at the University of Washington - Seattle since 2014. The original model was developed by Melissa Muradian as part of her M.S. under Dr. Trevor Branch. Dr. John Trochta made further modifications to the model, including incorporating an age-1 aerial school index and seroprevalence disease data during the course of his Ph.D. Joshua Zahner is currently utilizing BASA to develop a management strategy evaluation framework for PWS herring as a part of his M.S.

Cite as:

`Muradian ML, Branch TA, Moffitt SD, Hulson PJF (2017) Bayesian stock assessment of Pacific herring in Prince William Sound, Alaska. PLOS ONE 12(2): e0172153. https://doi.org/10.1371/journal.pone.0172153`

## Model Information
BASA is an age-structured population model using 10 distinct age classes: distinct classes for fish ages 0-8, and a plus-group class for fish ages 9+. Because herring recruit to the fishery at age 3, only age classes 3-9 are directly estimated by the model, while age classes 0-2 are back-calculated from the estimated number of age-3 fish. The model assumes a constant stock-recruitment relationship independent of annual spawning biomass, and a constant background natural mortality rate of M=0.25. Natural mortality is independently calculated in two half-seasons, summer and winter, and incorporates disease data. 

The model is reliant on data from several annual surveys performed in Prince William Sound. The mile-days-of-milt survey provides an estimate of total spawning output. The spawner age-composition survey provides estimates of age-3+ age composition. The hydroacoustic biomass surveys run by the Prince William Sound Science Center (PWSSC) provides an estimate of spawning biomass. And an age-1 school survey provides a rough estimation of the abundance of age-1 fish. The model also utilize three discontinued sources of data: an egg-deposition survey that provides an absolute index of spawning biomass, a second hydroacoustic survey (performed by Alaska Department of Fish and Game), and estimated age composition from the purse-seine fishery when it was active. Current work is attempting to integrate seroprevalence data from a disease survey into the model as well. 

## Dependencies

* ADMB-12.3

R Packages
* 'adnuts'   v. 1.1.2
* 'rstan'    v. 2.21.3
* 'snowfall' v. 1.84.6.1

## Running the Model
The actual stock assessment model runs in ADMB, though, in practice, is more freuqently run using the provided `run_basa.r` R script. `run_basa.R` acts as a wrapper that performs four important operations:

1. Calculates effective sample sizes (ESS) for each data source being provided to the model.
2. Generates initial values for MCMC sampler to begin at
3. Runs the stock assessment program using the NUTS sampling algorithm
4. Checks whether there were divergences between the sampling chains and whether all of the parameters were appropoiately estimated.

`run_basa.R` can be run directly from the command line or can be run line-by-line in RStudio or another IDE. 

Once the model has run succesfully, four R scripts are provided to visualize the results:
* `plotting/results_age_comp_predictions.R` generates a multi-panel bar plot displaying model estimated age compositions of the population for each year since 1980
* `plotting/results_model_survey_fits.R` generates a multi-panel plot displaying estaimte model fits to each of the five primary data sources.
* `plotting/results_recruitment_ssb.R` generates a two-panel plot displaying estimated recruitment to the fishery and estimated spawning biomass in each year since 1980.
* `plotting/results_output_for_management.R` generate a four panel plot summarising the model outputs for management use. The four plots include: annual recruitment, annual spawning biomass, annual fishing exploitation rate, and the posterior distribution of spawning biomass in the final year of the assessment. 

## Citations

The following scientific papers have made direct use of BASA in their work:

Monnahan, C. C., T. A. Branch, J. T. Thorson, I. J. Stewart, and C. S. Szuwalski. 2019. Overcoming long Bayesian run times in integrated fisheries stock assessments. ICES Journal of Marine Science 76:1477-1488.

Muradian, M. L., T. A. Branch, S. D. Moffitt, and P.J. F. Hulson. 2017. Bayesian stock assessment of Pacific herring in Prince William Sound, Alaska. PLOS ONE 12:e0172153.

Muradian, M. L., T. A. Branch, and A. E. Punt. 2019. A framework for assessing which sampling programmes provide the best trade-off between accuracy and cost of data in stock assessments. ICES Journal of Marine Science 76:2102-2113.

Trochta, J. T., T. A. Branch, A. O. Shelton, and D. E. Hay. 2020. The highs and lows of herring: A meta-analysis of patterns and factors in herring collapse and recovery. Fish and Fisheries 21:639-662.

Trochta, J. T., and T. A. Branch. 2021. Applying Bayesian model selection to determine ecological covariates for recruitment and natural mortality in stock assessment. ICES Journal of Marine Science 78:2875-2894.

Trochta, J. T., M. Groner, P. Hershberger, and T. A. Branch. 2022. A better way to account for disease in fisheries stock assessment: the powerful potential of seroprevalence data. Canadian Journal of Aquatic and Fisheries Sciences 79(4): 611-630



