# BASA model formulation


## Population dynamics

### Maturity-at-age ($\rho_{M,a}$)

Age 0-2 are assumed fully immature and age 5+ herring are assumed fully mature.
The maturity proportion of 3 and 4 herring are estimated parameters. Maturity of
age-3 herring ($\rho_{M,3} = \nu_3\nu_4$) is estimated as a fraction of the 
maturity proportion of age-4 herring ($\rho_{M,4} = \nu_4$).

$$ 
\rho_{M,a} = 
\begin{cases}
    0 & \text{ if } a < 3 \\
    \nu_3\nu_4 & \text{ if } a = 3 \\
    \nu_4 & \text{ if } a = 4 \\
    1 & \text{ if } a > 4
\end{cases}, \\
\text{where } \nu_3 \sim \text{Beta}(9, 11) \\ 
\text{ and } \nu_4 \sim \text{Beta}(18, 2) \\
$$

### Recruitment ($R_y$)

Recruitment ($R_y$) is modeled as a spawner-independent process wherein recruitment
deviates ($\epsilon_{Rec,y}$) vary log-normally about a bias-corrected mean 
recruitment parameter ($\bar{R}$).

$$
R_y = \bar{R} \ e^{\epsilon_{Rec,y} - 0.5\sigma_{Rec}^2}, \\
\text{where } \epsilon_{Rec,y} \sim \text{Normal}(0, \sigma_{Rec}^2) \\
\text{and } \sigma_{Rec} \sim \text{Uniform}(0.001, 2) \\
\text{and } \log(\bar{R}) \sim \text{Uniform}(2, 20)
$$

## Summary of estimated parameters

|          Parameter           |   Notation   |          Prior               |
| ---------------------------- | ------------ | ---------------------------- |
| Maturity proportion age-3 herring | $\nu_3: \rho_{M,3} = \nu_3\nu_4$ | $\text{Beta}(9,11)$ |
| Maturity proportion age-4 herring | $\nu_4: \rho_{M,4} = \nu_4$ | $\text{Beta}(18,2)$ |
| Mean recruitment (log-space) | $\log(\bar{R})$ |  $\text{Uniform}(2, 20)$ |
| Recruitment deviates | $\epsilon_{Rec,y}$ | $\text{Normal}(0, \sigma_{Rec}^2)$ |
| Recruitment deviates variance term | $\sigma_{Rec}$ | $\text{Uniform}(0.001, 2)$ |


## Other expressions in likelihood expressions

## Likelihoods


### L3

**Not scaled by number of years in data set**



### L5

**Not scaled by number of years in data set**

