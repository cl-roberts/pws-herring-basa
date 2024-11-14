## NUMBER OF SINGLE VALUE PARAMETERS
32

## —————————————————————————————————————————————————————————————————————————— ##
##                  DESIGN MATRIX FOR SINGLE VALUE PARAMETERS                 ##
##  Prior descriptions   Parameter values                                     ##
##  -0 uniform           (0,0)                                                ##
##  -1 normal            (p1=mu,p2=sig)                                       ##
##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
##  -3 beta              (p1=alpha,p2=beta)                                   ##
##  -4 gamma             (p1=alpha,p2=beta)                                   ##
## —————————————————————————————————————————————————————————————————————————— ##
##  init   lower   upper    est  prior                fun                     ##
## value   bound   bound    phz   type     p1    p2  type   # PARAMETER       ##
## —————————————————————————————————————————————————————————————————————————— ##
    0.08    0.00    5.00      3      0      0     0     0   # 1:  VHSV_age3_4_mort_93
    0.22    0.00    5.00      3      0      0     0     0   # 2:  ICH_age5_8_mort_93
    0.25    0.05    1.50     -2      0      0     0     0   # 3:  Z_0_8
    0.93    0.30    1.60      2      0      0     0     0   # 4:  Z_9
    0.60    0.01    0.90      3      3   9.00 11.00     0   # 5:  mat_par_1
    0.99    0.30    1.00      3      3  18.00  2.00     0   # 6:  mat_par_2
    0.60    0.01    0.90     -2      3   9.00 11.00     0   # 7:  mat_par_3
    0.90    0.30    1.00     -2      3  18.00  2.00     0   # 8:  mat_par_4
    3.66    3.00    5.00      4      0      0     0     0   # 9:  alpha_v
    2.83    1.00    7.00      4      0      0     0     0   # 10: beta_v
    4.00    3.00    5.00     -4      0      0     0     0   # 11: survey_vul_alpha
    2.40    1.00    7.00     -4      0      0     0     0   # 12: survey_vul_beta
    0.40    0.01    0.50     -2      0      0     0     0   # 13: egg_add
    5.87    2.30    9.00      2      0      0     0     0   # 14: logmdm_c
    0.33    0.01    0.90      5      1   0.33  0.10     0   # 15: m_add
   -0.38   -5.00    5.00      2      0      0     0     0   # 16: hydADFG_q
    0.30    0.01    0.70      5      1   0.30  0.08     0   # 17: hydADFG_add
   -0.21   -5.00    5.00      2      0      0     0     0   # 18: hydPWSSC_q
    0.32    0.01    0.60      5      1   0.32  0.08     0   # 19: hydPWSSC_add
    4.22   -5.00    8.00      1      0      0     0     0   # 20: log_q_juv
    1.00    0.01    4.00      5      0      0     0     0   # 21: overdisp_juv
    6.20    2.00   20.00     -1      0      0     0     0   # 22: log_MeanAge0
    0.00    0.00    2.00     -2      0      0     0     0   # 23: sigma_age0devs
    1.00    0.00    3.00     -5      0      0     0     0   # 24: sigma_mortdevs
    1.00   -3.00    5.00     -4      0      0     0     0   # 25: infec_vul_a50
    2.00   -2.00    8.00     -4      0      0     0     0   # 26: infec_vul_a95
    2.00   -3.00    5.00     -4      0      0     0     0   # 27: vhs_samp_vul_a50
    3.00   -2.00    8.00     -4      0      0     0     0   # 28: vhs_samp_vul_a95
    1.00   -3.00    5.00     -3      0      0     0     0   # 29: infec_vul_a50_ich
    2.00   -2.00    8.00     -3      0      0     0     0   # 30: infec_vul_a95_ich
    2.00    0.00    4.00     -3      0      0     0     0   # 31: ich_samp_vul_a50
    3.00    3.00    6.00     -3      0      0     0     0   # 32: ich_samp_vul_a95
## —————————————————————————————————————————————————————————————————————————— ##


## NUMBER OF VECTOR & MATRIX PARAMETERS 
13

## ROW IS A PARAMETER IN THE NEXT SECTION, AND COLUMN IS THE NUMBER OF ELEMENTS/COLS FOR EACH ROW BELOW EACH PARAMETER
1   2   7   5   # 1:  loginit_pop
1   4   7   1   # 2:  annual_age0devs
1   2   7   1   # 3:  beta_age0
1   2   7   1   # 4:  beta_mortality
1   4   7   1   # 5:  annual_mortdevs
1   2   7   1   # 6:  beta_age0_offset
1   2   7   1   # 7:  beta_mortality_offset
1   2   7   1   # 8:  sigma_age0covar
1   2   7   1   # 9:  sigma_morcovar
1   2   7   1   # 10: inf_prob_vhs
1   2   7   1   # 11: rec_prob_vhs
1   2   7   1   # 12: inf_prob_ich
1   2   7   1   # 13: rec_prob_ich

## —————————————————————————————————————————————————————————————————————————— ##
##                CONTROLS FOR VECTOR & MATRIX PARAMETERS                     ##
## —————————————————————————————————————————————————————————————————————————— ##

## 1: loginit_pop
   1                                    # Dimensions (1=Vector, 2=Matrix)
   1    5                               # Start & end of each dimension for initial value(s) listed below
   3.00    8.00   1   0   0   0   0     # lower bound, upper bound, phz, prior type, p1, p2, fun type
   6.35  5.66  5.92  6.74  4.74         # Initial value(s)

## 2: annual_age0devs
   2                                    
   1   1   1   1                        # While only a single value, indices are set internally (within TPL)
   -10.00 10.00   1   0   0   0   0     # These have to be really wide for the MSE, or else the hessian isn't P&D
   0

## 3: beta_age0
   1
   1   1
   -10.00 10.00 -1   0   0   0   0
   0

## 4: beta_mortality
   1
   1   1
   -30.00 30.00  2   0   0   0   0
   0.2

## 5: annual_mortdevs
   2
   1   1   1   1
   -20.00 20.00  -2   0   0   0   0
   0

## 6: beta_age0_offset
   1
   1   1
   -5.00  -5.00  -1   0   0   0   0
   0

## 7: beta_mortality_offset
   1
   1   1
   -5.00  -5.00  -1   0   0   0   0
   0

## 8: sigma_age0covar
   1
   1   1
   0.01    2.00  -5   0   0   0   0
   0.5

## 9: sigma_morcovar
   1
   1   1
   0.01    2.00  -5   0   0   0   0
   0.5

## 10: inf_prob_vhs
   1
   29   29
#  0.00    1.00   -3   3   1.5  5   0  # Use commented out row for estimating these parameters
   0.00    1.00   -3   0   0   0   0
   0.0

## 11: rec_prob_vhs
   1
   29   29
#  0.0    0.99   -3   3   4   6   0
   0.0    0.99   -3   0   0   0   0
   0.0

## 12: inf_prob_ich
   1
   25   25
#  0.0    1.00   -3   3   2   8   0
   0.0    1.00   -3   0   0   0   0
   0.0

## 13: rec_prob_ich
   1
   25   25
#  0.0    0.99   -3   3   5   5   0
   0.0    0.99   -3   0   0   0   0
   0.0

## —————————————————————————————————————————————————————————————————————————— ##
