# Appendix B.-- Prince William Sound herring spring pre-fishery mature biomass forecast methods for Bayesian age-structured assessment (BASA)

<!-- revisit projecting N_y_a foreward rather than using 10-year median -->

> [!note]
> Notation: $(*)$ denotes forecast quantity.

The Bayesian age-structured assessment (BASA) model is used to forecast the spring 
pre-fishery mature biomass of Prince William Sound herring immediately prior to the 
onset of spring fisheries. This quantity is also sometimes referred to as pre-fishery 
spawning biomass. In summary, a posterior distribution for spring pre-fishery mature biomass 
is estimated by projecting the number of age 2--9+ fish in fall using winter mortality 
based on a ten-year average and accounting for fall/winter catch. The numbers 
of age-2 fish are forecasted using a ten-year average of model estimates and the 
numbers of age 3--9+ fish are estimated by the model. Numbers-at-age are converted 
to mature biomass using ten-year averages for spawning weight-at-age and model-estimated 
maturity.

## Numbers-at-age forecast

We forecast age-3 spring recruitment based on a geometric mean of the
previous 10 recruitment sizes estimated by the model. Let year $Y$ denote the final year 
of the model time series and, thus, year $Y+1$ is the forecast year. During the model fitting 
process, we develop a posterior distribution for the age-3 recruitment forecast
by generating a sample of age-2 fish for year $Y$, applying summer mortality in 
year $Y$, subtracting the age-2 fish caught in the year $Y$ food/bait fishery, 
and applying year-$Y$ winter mortality. The number of age-2 fish in year $Y$ is
sampled from a lognormal distribution parameterized by the mean and standard
deviation of age-2 fish, in log-space, estimated for the previous ten years:

$$
r^*_{\text{age} 2} \sim \text{lognormal} \left( \mu_{\text{age} 2}, \sigma_{\text{age} 2}\right) 
$$

where 

$$
\mu_{\text{age} 2} = \text{mean}\{\log N_{y, 2}\} \ \ \ \text{and} \ \ \  
\sigma_{\text{age} 2} = \text{SD}\{\log N_{y, 2}\} \ \ \
\forall \ y \ \in \ [Y-9, Y],
$$

and thus the forecasted numbers-at-age matrix in year $Y$ is

$$
N^*_{Y+1, a} = 
  \begin{cases} 
    \left(r^*_{\text{age} 2} S^{1*}_{Y, a-1} - C^4_{Y, a-1}\right) S^{2*}_{Y, a-1} & \text{if } a = 3 
    \\ \\
    \left((N_{Y, a-1}-\hat{C}^S_{Y, a-1}) S^{1*}_{Y, a-1} - C^4_{Y, a-1}\right) S^{2*}_{Y, a-1} & \text{if } a \geq 4 \text{ and } a \leq 8 
    \\ \\
    \left((N_{Y, a-1}-\hat{C}^S_{Y, a-1}) S^{1*}_{Y, a-1} - C^4_{Y, a-1}\right) S^{2*}_{Y, a-1}  
    + \left((N_{Y, a}-\hat{C}^S_{Y, a}) S^{1*}_{Y, a} - C^4_{Y, a}\right) S^{2*}_{Y, a}  & \text{if } a = 9+.
  \end{cases}
$$

Note that $y$ indexes year, $a$ indexes age, and $(*)$ denotes a forecasted quantity. 
The summer ($S^{1*}_{Y, a}$) and winter ($S^{1*}_{Y, a}$) survival forecasts are based 
on 10-year averages. Catches from the most recent spring and fall commercial fisheries 
are given by $\hat{C}^S_{Y, a}$ and $C^4_{Y, a}$, respectively.


## Spring pre-fishery mature biomass forecast 

The spring pre-fishery mature biomass forecast is given by 

$$
\tilde{B}^*_{Y+1} = \sum_{a = 0}^{9+} \rho_a \cdot \bar{W}_{Y+1, a} \cdot N^*_{Y+1, a}
$$

where the weight-at-age used in the forecast ($\bar{W}_{Y+1, a}$) is based on a ten-year
average of sampled weights of spawning fish and maturity ($\rho_a$) is a model-estimated 
parameter.
