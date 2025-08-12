# BASA model formulation

## Data sets

Let $Y$ denote the complete set of years fit by the model, beginning in 1980 and ending 
after $n_Y$ years ($Y = [1980, 1980+n_Y]$). $A$ is the set of age classes where
age-9 is the plus group ($A = [0, 9+]$). The datasets fit by BASA are:


| Data type | Units | Symbol | Years |
| :----- | :----: | :---: | :---: |
| Gillnet catch-at-age | millions | $C^2_{y, a}$ | $[1980, 1998]$ |
| Pound utilization catch-at-age | millions | $C^3_{y, a}$ | $[1980, 1999]$ |
| Food/bait catch-at-age | millions | $C^4_{y, a}$ | $[1980, 1998]$ |
| Fecundity-at-age | no. of eggs per female | $f_{y, a}$ | $[1984, 1993]$ |
| Weight-at-age of spawning herring | $\frac{\text{mt}}{\text{million fish}}$ | $w_{y, a}$ | $[1980, 1980+n_Y]$ |
| Purse-seine age-composition | proportion | $\Theta^1_{y, a}$ | $[1980, 1998]$ |
| Spawner survey age-composition | proportion | $\Theta^\text{Sp}_{y, a}$ | $[1980, 1980+n_Y]$ |
| Female spawners | proportion | $\rho^f_y$ | $[1980, 1980+n_Y]$ |
| Total annual purse-seine yield | mt | $\Omega^1_{y}$ | $[1980, 1998]$ |
| Eggs deposited | trillions | $E_y$ | $[1984, 1997]$ |
| C.V. for eggs deposited | -- | $\sigma^E_y$ | $[1984, 1997]$ |   
| ADFG hydroacoustic survey biomass | mt | $H^A_y$ | $[2005, 2009]$ |   
| PWSSC hydroacoustic survey biomass | mt | $H^P_y$ | $[1993, 2020]$ |   
| C.V. for PWSSC hydroacoustic biomass | -- | $\sigma^{H^P}_y$ | $[1993, 2020]$ |   
| Milt | $\frac{\text{miles}}{\text{day}}$ | $T_y$ | $[1980, 1980+n_Y]$ |
| Juvenile schools | no. age-1 schools (small school equivalents) | $J_y$ | $[2010, 1980+n_Y]$ |
| *Ich.* prevalence (pre-2007 survey) | Proportion | $D^{Ich. 1}_y$ | $[1994, 2006]$ |
| *Ich.* prevalence (post-2007 survey)[^1] | Proportion | $D^{Ich. 2}_y$ | $[2007, 1980+n_Y]$ |
| VHS prevalence | Proportion | $D^{VHS}_y$ | $[1994, 1980+n_Y]$ |

[^1]: *Ichthyophonus* prevalence sampling (?) methods changed beginning in 2007


## Model formulation

### Population model

Let $b$ index half-year, where $b=1$ denotes summer and $b=2$ winter. The population
dynamics modeled by BASA are:

| Description | Equation | 
| :----- | :-----: |
*Catch, millions of fish*
| Estimated total purse-seine catch | $\hat{C}^1_y = \frac{\Omega^1_{y}}{\sum_{a \in A} (\hat{\Theta}^1_{y, a} w_{y, a})}$ |
| Spring removals | $\hat{C}^S_{y, a} = \hat{\Theta}^1_{y, a} \hat{C}^1_y + C^2_{y, a} + \rho_k C^3_{y, a}$ |
| Spawn removals | $\tilde{C}^S_{y, a} = \hat{\Theta}^1_{y, a} \hat{C}^1_y + C^2_{y, a} + C^3_{y, a}$ |
*Half-year survival, proportion*
| Ages 0+, 1980 - 1991 | $S^b_{y, a} = e^{-0.5 \bar{M}}$ |
| Ages 0-2, 1992 - 1980+n_Y | $S^b_{y, a} = e^{-0.5 \bar{M}}$ |
| Ages 3-4, 1992 - 1993 | $S^b_{y, a} = e^{-0.5 \bar{M} - M^{VHS}_{1992-1993}}$ |
| Ages 3-4, 1994 - 1980+n_Y | $S^b_{y, a} = e^{-0.5 \bar{M} - \beta_3 D^{VHS}_y}$ |
| Ages 5-8, 1992 - 1993 | $S^b_{y, a} = e^{-0.5 \bar{M} - M^{Ich.}_{1992-1993}}$ |
| Ages 5-8, 1994 - 2006 | $S^b_{y, a} = e^{-0.5 \bar{M} - \beta_1 D^{Ich. 1}_y}$ |
| Ages 5-8, 2007 - 1980+n_Y | $S^b_{y, a} = e^{-0.5 \bar{M} - \beta_2 D^{Ich. 2}_y}$ |
| Ages 9+, 1980 | $S^b_{1980, 9+} = e^{-0.5 \bar{M}_{9+}}$ |
| Ages 9+, 1981 - 1980+n_Y | $S^b_{y, 9+} = S^b_{y-1, 9+} \left( \frac{S^b_{y, 8}}{S^b_{y-1, 8}} \right)$ |
*Recruitment, millions*
| Annual recruitment (age-0) | $R_y = \bar{R} e^{\epsilon^{rec}_y}$ |
*Selectivity, proportion*
| Purse-seine age-specific selectivity | $V_a = \frac{1}{1 + e^{-\beta^V (a-\alpha^V)}}$ |
*Numbers-at-age, millions*
| Pre-fishery total abundance, age 0-2 | $N_{y+1, a+1} = N_{y, a} S^1_{y, a} S^2_{y, a}$ |
| Pre-fishery total abundance, ages 3-8 | $N_{y+1, a+1} = [(N_{y, a}-\hat{C}^S_{y, a}) S^1_{y, a} - C^4_{y, a}] S^2_{y, a}$ |
| Pre-fishery total abundance, age 9+ | $N_{y+1, 9+} = [(N_{y, 8}-\hat{C}^S_{y, 8}) S^1_{y, 8} - C^4_{y, 8}] S^2_{y, 8} + [(N_{y, 9+}-\hat{C}^S_{y, 9+}) S^1_{y, 9+} - C^4_{y, 9+}] S^2_{y, 9+}$ |
| Post-fishery spawning abundance | $\tilde{N}_{y, a} = \rho^M_a (N_{y, a} - \tilde{C}^S_{y, a})$ |
*Biomass, metric tons*
| Pre-fishery total biomass | $B_y = \sum_{a \in A} N_{y, a} W_{y, a}$ |
| Pre-fishery mature biomass | $\tilde{B}_y = \sum_{a \in A} \rho^M_a N_{y, a} W_{y, a}$ |
| Post-fishery spawning biomass | $\tilde{B}^{post}_y = \sum_{a \in A} \tilde{N}_{y, a} W_{y, a}$ |


### Data model

Other expressions used in the likelihood components of BASA are:

| Description | Equation | 
| :----- | :-----: |
| Estimated purse-seine age composition, proportion | $\hat{\Theta}^1_{y, a} = \frac{V_a N_{y, a}}{\sum_{a \in A} (V_a N_{y, a})}$ |
| Estimated spawner survey age composition, proportion | $\hat{Sp}^1_{y, a} = \frac{\rho^M_a N_{y, a}}{\sum_{a \in A} (\rho^M_a N_{y, a})}$ |
| Estimated egg deposition, trillions | $\hat{E}_y = 10^{-6} \rho^f_y \sum_{a \in A} (f_{y, a} \tilde{N}_{y, a})$ |
| Estimated ADFG hydroacoustic biomass, metric tons | $\hat{H}^A_y = \tilde{B}_y \cdot e^{q^{H^A}}$ |
| Estimated PWSSC hydroacoustic biomass, metric tons | $\hat{H}^P_y = \tilde{B}_y \cdot e^{q^{H^P}}$ |
| Estimated milt, mile-days | $\hat{T}_y = \frac{(1-\rho^f_y) \tilde{B}^{post}_y}{e^{q^T}}$ |
| Estimated juvenile schools, small-school equivalent | $\hat{J}_y = N_{y, 1} \cdot e^{q^J}$ |

## Model parameters 

The total number of parameters estimated by BASA are:

| Parameter | Symbol | Prior |
| :----- | :---: | :---: |
*Fixed parameters*
| Background natural mortality, ages 0-8 | $\bar{M}$ | -- |
| Egg deposition additional error | $\sigma^{E+}$ | -- |
| Mean recruitment | $\bar{R}$ | -- |
| Pound mortality | $\rho_k$ | -- |
*Initial population size*
| Age-1 abundance in 1980, log-link | $\eta_{1980, 1}; N_{1980, 1} = e^{\eta_{1980, 1}}$ | uniform |
| Age-2 abundance in 1980, log-link | $\eta_{1980, 2}; N_{1980, 2} = e^{\eta_{1980, 2}}$ | uniform |
| Age-3 abundance in 1980, log-link | $\eta_{1980, 3}; N_{1980, 3} = e^{\eta_{1980, 3}}$ | uniform |
| Age-4 abundance in 1980, log-link | $\eta_{1980, 4}; N_{1980, 4} = e^{\eta_{1980, 4}}$ | uniform |
| Age-5 abundance in 1980, log-link | $\eta_{1980, 5}; N_{1980, 5} = e^{\eta_{1980, 5}}$ | uniform |
*Recruitment deviates*
| Annual recruitment deviations, log-link | $\epsilon^{rec}_y; N_{y, 0} = \bar{R} e^{\epsilon^{rec}_y}$ | uniform |
*Mortality parameters*
| Background natural mortality, ages 9+ | $\bar{M}_{9+}$ | uniform |
| Additional ages 3-4 mortality, 1992-1993 | $\bar{M}^{VHS}_{1992-1993}$ | uniform |
| Additional ages 5-8 mortality, 1992-1993 | $\bar{M}^{Ich.}_{1992-1993}$ | uniform |
| Additional mortality due to VHS, 1994 - 1980+n_Y | $\beta^1$ | uniform |
| Additional mortality due to *Ich.*, 1994 - 2006 | $\beta^2$ | uniform |
| Additional mortality due to *Ich.*, 2007 - 1980+n_Y | $\beta^3$ | uniform |
*Selectivity parameters*
| Purse-seine gear selectivity | $\alpha^V$ | uniform |
| Purse-seine gear selectivity | $\beta^V$ | uniform |
*Maturity parameters*
| Age-3 maturity | $\nu_3; \rho^M_3 = \nu_3 \rho^M_4$ | beta |
| Age-4 maturity | $\rho^M_4$ | beta |
*Data model parameters*
| ADFG Hydroacoustic catchability coefficient | $q^{H^A}$ | uniform |
| ADFG Hydroacoustic biomass CV | $\sigma^{H^A}$ | normal |
| PWSSC Hydroacoustic catchability coefficient | $q^{H^P}$ | uniform |
| PWSSC Hydroacoustic biomass additonal error | $\sigma^{H^P+}$ | normal |
| Milt catchability coefficient | $q^T$ | uniform |
| Milt CV | $\sigma^T$ | normal |
| Juvenile schools catchability coefficient | $q^J$ | uniform |
| Juvenile schools overdispersion | $\phi^J$ | uniform |


## Likelihood expressions

### Effective sample sizes

Due to herring schooling behavior and gear selectivity, effective sample sizes 
were estimated to weight the likelihood expressions for purse-seine age composition
($L_1$) and the ADFG spawner survey age-composition ($L_2$) datasets. A iterative 
reweighting procedure is conducted for estimating effective sample sizes wherein 
BASA is initially using raw sample sizes for the age-composition datasets, and
the effective sample size is estimated using

$$
\tilde{Z}^1_y = \sum_{a=3}^{9+} \frac{\hat{\Theta}^1_{y, a} (1 - \hat{\Theta}^1_{y, a})}{(\Theta^1_{y, a} - \hat{\Theta}^1_{y, a})^2}

\hspace{.5cm}
\text{and}
\hspace{.5cm}

\tilde{Z}^{Sp}_y = \sum_{a=3}^{9+} \frac{\hat{\Theta}^{Sp}_{y, a} (1 - \hat{\Theta}^{Sp}_{y, a})}{(\Theta^{Sp}_{y, a} - \hat{\Theta}^{Sp}_{y, a})^2}.
$$

These effective sample sizes are then used to re-fit the model, and the procedure
repeats until the harmonic mean (across years) of the ratio of the effective to raw 
sample size converges.

### Likelihood components

Let $n_L$ denote the number of likelihood components. $Y_i$ gives the set of 
years for the dataset in likelihood $i$ where $i \in [1, n_L]$. The cardinality 
of $Y_i$ is given by $n_i$. The negative log-likelihoods used in the BASA objective 
function are:

| Likelihood component | Assumed probability distribution | Form |
| :----- | :---: | :-----: |
| Complete expression | -- | $L = \sum_{i=1}^{n_L} L_i$ |
| Purse-seine age-composition | multinomial | $L_1 = -\sum_{y \in Y_1} \left[ \tilde{Z}^1_y \sum_{a \in A} \Theta^1_{y, a} \log\left( \frac{\hat{\Theta}^1_{y, a}}{\Theta^1_{y, a}} \right) \right]$ |
| Spawner survey age-composition | multinomial | $L_2 = -\sum_{y \in Y_2} \left[ \tilde{Z}^{Sp}_y \sum_{a \in A} \Theta^{Sp}_{y, a} \log\left( \frac{\hat{\Theta}^{Sp}_{y, a}}{\Theta^{Sp}_{y, a}} \right) \right]$ |
| Egg deposition | lognormal | $L_3 = \sum_{y \in Y_3} \left[ \log(\sigma^{L_3}_y) + \frac{\left(\log(\hat{E}_y) - \log(E_y)\right)^2}{2 (\sigma^{L_3}_y)^2} \right]$ |
| Total variance for $L_3$ | -- | $\left(\sigma^{L_3}_y\right)^2 = \left(\sigma^E_y\right)^2 + \left(\sigma^{E+}\right)^2$ |
| ADFG hydroacoustic biomass | lognormal | $L_4 = n_4 \log(\sigma^{H^A}) + \frac{1}{(\sigma^{H^A})^2} \sum_{y \in Y_4} \left[ \left(\log(\hat{H}^A_y) - \log(H^A_y)\right)^2 \right]$ |
| PWSSC hydroacoustic biomass | lognormal | $L_5 = \sum_{y \in Y_5} \left[ \log(\sigma^{L_5}_y) + \frac{\left(\log(\hat{H}^p_y) - \log(H^p_y)\right)^2}{2 (\sigma^{L_5}_y)^2} \right]$ |
| Total variance for $L_5$ | -- | $\left(\sigma^{L_5}_y\right)^2 = \left(\sigma^{H^P}_y\right)^2 + \left(\sigma^{H^P+}\right)^2$ |
| Milt | lognormal | $L_6 = n_6 \log(\sigma^T) + \frac{1}{(\sigma^T)^2} \sum_{y \in Y_6} \left[ \left(\log(\hat{T}_y) - \log(T_y)\right)^2 \right]$ |
| Juvenile schools [^2] | negative binomial | $L_7 = - \sum_{y \in Y_7} \left[ \log \Gamma\left(J_y+m_y\right) + \log \Gamma\left(m_y\right) + \log \Gamma\left(1+J_y\right) - m_y \log(p_y) - J_y \log(1-p_y) \right]$ |
| Total variance for $L_7$ | -- | $\left(\sigma^{L_7}_y\right)^2 = \hat{J}_y + \frac{\hat{J}_y^2}{\phi^J}$ |

[^2]: $\Gamma$ denotes the gamma function; $p_y = \frac{\hat{J}_y}{\left(\sigma^{L_7}_y\right)^2}$ and $m_y = \hat{J}_y \cdot \frac{p_y}{1-p_y}$