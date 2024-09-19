// ========================================================================== //
//                                                                            //
//       Bayesian age-structured model for Prince William Sound herring       //
//                                                                            //
//                 Written in Template Model Builder (TMB)                    //
//                                                                            //
//                            UNDER DEVELOPMENT                               //
//                                Aug  2023                                   //
//                                                                            //
//                                 AUTHORS                                    //
//                                CL Roberts                                  //
//                           cl.roberts@alaska.gov                            //
//                                                                            //
//               Adapted from BASA model developed in ADMB at UW              //
//              by John Trochta, Trevor Branch, and Joshua Zahner.            //
//                                                                            //
//            The original ADMB model was built on code developed by          //
//             Melissa Muradian (UW) and adapted from an Excel-based          //
//                      model by Steven Moffitt (ADF&G)                       //
//                                                                            //
// -------------------------------------------------------------------------- //
//                                                                            //
// Program file:  PWS_ASA_tmb.cpp                                             //
// Input Data files                                                           //
//        PWS_ASA.dat:             Model input (surveys, catches, etc.)       //
//        PWS_ASA(phases):         Parameter phases                           //
//        PWS_ASA(ESS).ctl:        Effective sample sizes iteratively         //
//                                 estimated (in R) and used in obj function  //
//        pws_asa.PIN:             Included to test different starting values //
//                                                                            //
// Output files:                                                              //
//        PWS_ASA.rep: Customizable file in the REPORT_SECTION, useful for    //
//                     trouble-shooting                                       //
//        PWS_ASA.std: Default ADMB report file with Hessian derived SE's     //
//        PWS_ASA.par: ADMB parameter estimates                               //
//                                                                            //
// -------------------------------------------------------------------------- //
//                                                                            //
// Refer to model documentation                                               //
//                                                                            //
// ========================================================================== //


// -------------------------------------------------------------------------- //
//                             Globals Section                                //
// -------------------------------------------------------------------------- //

// ---- load libraries ---- //

#include <TMB.hpp>


// ---- utility functions ---- //

// density function for uniform distribution
template <class Type>
Type dunif(Type x, Type lwr, Type upr, int give_log=0)
{
    Type dens = 0;

    if (lwr >= upr) {
        throw std::invalid_argument("upper bound is not larger than lower bound");
    }

    if (x < lwr || upr < x) {
        dens = 0;
    } else if (lwr <= x && x <= upr) {
        dens = 1 / (upr-lwr);
    } 

    if (give_log) return log(dens); else return dens;
}

template <class Type>
Type posfun(Type x, Type eps, Type pen)
{
    if (x>=eps) {
        return x;
    } else {
        pen += .01*pow(x-eps, 2);
        return eps / (2-x/eps);
    }
}


// ---- define objective function ---- //

template<class Type>
Type objective_function<Type>::operator() ()
{

    // ------------------------------------------------------------------ //
    //                          Data Section                              //
    // ------------------------------------------------------------------ //

    // ---- model/PWS_ASA(ESS).ctl ---- //
    DATA_IVECTOR(seine_ess);                // purse-seine effective sample sizes
    DATA_IVECTOR(spawn_ess);                // spawn survey effective sample sizes
    DATA_IVECTOR(vhsv_ess);                 // VHSV effective sample sizes

    // ---- model/PWS_ASA.dat ---- //
    DATA_INTEGER(nyr);                      // length of model data time series
    // DATA_INTEGER(nyr_tobefit);              // number of years fit by the model
    DATA_INTEGER(nage);                     // number of age classes
    DATA_MATRIX(waa);                        // weight-at-age matrix
    DATA_MATRIX(fecundity);                  // fecundity-at-age matrix
    DATA_MATRIX(pound_catch);                // pound catch-at-age
    DATA_MATRIX(foodbait_catch);             // food/bait catch-at-age matrix
    DATA_MATRIX(gillnet_catch);              // gillnet catch-at-age matrix
    DATA_VECTOR(seine_yield);               // total purse-seine yield over time
    DATA_VECTOR(perc_female);               // female proportion of spawners 
    DATA_VECTOR(mdm);                       // mile-days milt survey index
    DATA_VECTOR(egg);                       // egg deposition survey index
    DATA_VECTOR(egg_se);                    // egg deposition survey se's
    // DATA_INTEGER(adfg_hydro_year_start);    // index first year of ADFG hydro survey
    DATA_IVECTOR(adfg_hydro);               // ADFG hydroacoustic survey index    
    // DATA_INTEGER(pwssc_hydro_year_start);   // index first year of PWSSC hydro survey
    DATA_IVECTOR(pwssc_hydro);              // PWSSC hydroacoustic survey index
    DATA_VECTOR(pwssc_hydro_se);            // PWSSC hydroacoustic survey se's
    DATA_MATRIX(seine_age_comp);             // purse-seine age compositions
    DATA_MATRIX(spawn_age_comp);             // spawn survey age compositions
    DATA_IVECTOR(juvenile_survey);          // aerial juvenile school survey index

    // ---- model/PWS_ASA(covariate).ctl ---- //
    // mortality covariates:
        // | I. hoferi pre-2007  |  I. hoferi post-2007  |  VHSV |
    DATA_MATRIX(disease_covs);
    DATA_IMATRIX(mort_age_impact);           // which age does the mort cov impact?

    // ---- model/agecomp_samp_sizes.txt ---- //
    DATA_IVECTOR(seine_sample_size);        // purse seine raw sample sizes
    DATA_IVECTOR(spawn_sample_size);        // spawn survey raw sample sizes
    DATA_IVECTOR(vhsv_sample_size);         // antibody raw sample sizes
    DATA_IVECTOR(ich_sample_size);          // I. hoferi raw sample sizes


    // ------------------------------------------------------------------ //
    //                        Parameter Section                           //
    // ------------------------------------------------------------------ //

    // ---- fixed parameters ---- //
    PARAMETER(pk);                      // pound spawn-on-kelp mortality proportion
    PARAMETER(egg_add);                 // egg deposition additional variance term
    PARAMETER(Z_0_8);                   // background mortality, ages 0-8
    PARAMETER(log_MeanAge0);            // mean recruits, log-space
    PARAMETER(sigma_age0devs);          // variance term for recruitment deviates

    // ---- estimated parameters ---- //
    PARAMETER_VECTOR(annual_age0devs);  // recruitment deviates
    PARAMETER(log_juvenile_q);          // aerial juvenile survey scalar, log-space
    PARAMETER_VECTOR(loginit_pop);      // 1980 numbers-at-age, ages 1-5, log-space

    PARAMETER(Z_9);                     // background mortality, plus group (9+)
    PARAMETER_VECTOR(beta_mortality);   // disease covariate effect on survival
    PARAMETER(logmdm_c);                // milt scalar, log-space                        
    PARAMETER(adfg_hydro_q);            // hydroacoustic scalar (ADFG)
    PARAMETER(pwssc_hydro_q);           // hydroacoustic scalar (PWSSC)

    PARAMETER(VHSV_age3_4_mort_93);     // additional mortality 1992-93 ages 3-4 (VHSV)
    PARAMETER(ICH_age5_8_mort_93);      // additional mortality 1992-93 ages 5-8 (I. hoferi)
    PARAMETER(mat_age3);                // maturity prop age 3 (as frac of age 4 prop)
    PARAMETER(mat_age4);                // maturity proportion age 4

    PARAMETER(seine_selex_alpha);       // purse-seine selectivity logistic alpha parameter
    PARAMETER(seine_selex_beta);        // purse-seine selectivity logistic beta parameter

    PARAMETER(milt_add_var);            // milt additional variance
    PARAMETER(adfg_hydro_add_var);      // hydroacoustic additional variance (ADFG)
    PARAMETER(pwssc_hydro_add_var);     // hydroacoustic additional variance (PWSSC)
    PARAMETER(juvenile_overdispersion); // aerial juvenile additional variance


    // ------------------------------------------------------------------ //
    //                        Procedure Section                           //
    // ------------------------------------------------------------------ //

    // ---- convert some parameters from log-space ---- //
    Type MeanAge0 = exp(log_MeanAge0);                  // mean recruitment
    vector<Type> init_pop = exp(loginit_pop);           // initial population vector

    // --------------------- POPULATION DYNAMICS ------------------------ //

    // ---- survival from natural mortality ---- //
    // calculate half-year survival using background mortality, plus-group
    // background mortality, excess mortality in 1992-1993 for the different
    // age groups, and 3 disease covariates

    // 2 matrices (summer and winter) will be used to store half-year survival 
    matrix<Type> summer_survival(nyr, nage);
    summer_survival = summer_survival.setZero();
    matrix<Type> summer_mortality_effect(nyr, nage);
    summer_mortality_effect = summer_mortality_effect.setZero(); 
    matrix<Type> winter_survival(nyr, nage);
    winter_survival = winter_survival.setZero();
    matrix<Type> winter_mortality_effect(nyr, nage);
    winter_mortality_effect = winter_mortality_effect.setZero(); 
    int index_1992 = 1992-1980;                 // index for 1992
    int n_covs = beta_mortality.size();

    matrix<Type> disease_covs_calc(nyr, n_covs);

    for (int y = 0; y < nyr; y++) {
        for(int a = 0; a < nage; a++) {
            for(int b = 0; b < n_covs; b++){
                if (disease_covs(y, b) == -9) {
                    disease_covs_calc(y, b) = 0;
                } else if (disease_covs(y, b) != -9) {
                    disease_covs_calc(y, b) = disease_covs(y, b);
                }
                if(y == nyr-1){
                    summer_mortality_effect(y, a) += beta_mortality(b)*disease_covs_calc(y, b)*mort_age_impact(a, b);
                } else {
                    summer_mortality_effect(y, a) += beta_mortality(b)*disease_covs_calc(y, b)*mort_age_impact(a, b);
                    winter_mortality_effect(y+1, a) += beta_mortality(b)*disease_covs_calc(y, b)*mort_age_impact(a, b);
                }             
            }
        }
    }

    for (int y = 0; y < nyr; y++) {
        for (int a = 0; a < (nage-1); a++) {
            // ages 0-8
            summer_survival(y, a) = exp(-0.5*Z_0_8);
            winter_survival(y, a) = exp(-0.5*Z_0_8);
            for(int b = 0; b < n_covs; b++) {
                // if (disease_covs(y, b) == -9) continue; 
                summer_survival(y, a) *= exp(-summer_mortality_effect(y, a));
                winter_survival(y, a) *= exp(-winter_mortality_effect(y, a));
            }
            if (a >= 3 && y == index_1992) {
                if (a <= 4) {
                    summer_survival(y, a) *= exp(-VHSV_age3_4_mort_93);
                } else if (a > 4) {
                    summer_survival(y, a) *= exp(-ICH_age5_8_mort_93);
                }
            } else if (a >= 3 && y == index_1992+1) {
                if (a <= 4) {
                    summer_survival(y, a) *= exp(-VHSV_age3_4_mort_93);
                    winter_survival(y, a) *= exp(-VHSV_age3_4_mort_93);
                } else if (a > 4) {
                    summer_survival(y, a) *= exp(-ICH_age5_8_mort_93);
                    winter_survival(y, a) *= exp(-ICH_age5_8_mort_93);
                }              
            } else if (a >= 3 && y == index_1992+2) {
                if (a <= 4) {
                    winter_survival(y, a) *= exp(-VHSV_age3_4_mort_93);
                } else if (a > 4) {
                    winter_survival(y, a) *= exp(-ICH_age5_8_mort_93);
                }    
            }
        } 
        // plus group (9+)
        int a = nage-1;
        if (y == 0) {
            summer_survival(y, a) = exp(-0.5*Z_9);
            winter_survival(y, a) = exp(-0.5*Z_9);
            for(int b = 0; b < n_covs; b++) {
                // if (disease_covs(y, b) == -9) continue; 
                summer_survival(y, a) *= exp(-summer_mortality_effect(y, a));
                winter_survival(y, a) *= exp(-winter_mortality_effect(y, a));
            }                             
        } else if (y > 0) {
            summer_survival(y, a) = summer_survival(y-1, a) * (summer_survival(y, a-1)/summer_survival(y-1, a-1));
            winter_survival(y, a) = winter_survival(y-1, a) * (winter_survival(y, a-1)/winter_survival(y-1, a-1));
        }
    }

    // penalize survival higher than 1

    Type max_surv = .99;
    Type winter_surv_pen = 0.0;
    Type summer_surv_pen = 0.0;

    for (int y = 0; y < nyr; y++) { 
        for (int a = 0; a < nage; a++) {
            if (winter_survival(y, a) > max_surv) {
                winter_surv_pen += .01*pow(winter_survival(y, a) - max_surv, 2);
                winter_survival(y, a) = max_surv;
            }
            if (summer_survival(y, a) > max_surv) {
                summer_surv_pen += .01*pow(summer_survival(y, a) - max_surv, 2);
                summer_survival(y, a) = max_surv;
            }
        }
    }


    // ---- calculate seine selectivity ---- //
    vector<Type> seine_selex(nage);                     // seine selectivity
    for (int a = 0; a < nage; a++) {
        if (a < 3) {
            seine_selex(a) = 0;
        } else if (a >= 3) {
            seine_selex(a) = 1 / (1 + exp(-seine_selex_beta*(a - seine_selex_alpha)));
        }
    }

    // ---- maturation schedule ---- //
    vector<Type> maturity(nage);                        // age-varying maturity proportion
    for (int a = 0; a < nage; a++) {
        if (a < 3) {
            maturity(a) = 0;
        } else if (a == 3) {
            maturity(a) = mat_age3 * mat_age4;
        } else if (a == 4) {
            maturity(a) = mat_age4;
        } else if (a > 4) {
            maturity(a) = 1;
        }
    }

    // ---- recruitment ---- //
    vector<Type> age_0(nyr);
    for (int y = 0; y < nyr-1; y++) {
        age_0(y) = MeanAge0 * exp(annual_age0devs(y) - 0.5*pow(sigma_age0devs, 2));
        // age_0(y) = exp(annual_age0devs(y));
    }
    // note: the TPL model sets the final recruitment deviation to 0, effectively
    // forecasting the age-0 herring in the terminal year of the model to be
    // the MeanAge0 (currently a fixed parameter) 
    age_0(nyr-1) = MeanAge0;                                       


    // ---- initialize pop dynamics objects and fill in first year ---- //

    // pre-fishery numbers-at-age matrix, first row
    matrix<Type> N_y_a(nyr, nage);
    N_y_a = N_y_a.setZero(nyr, nage);              
    N_y_a(0, 0) = age_0(0);
    for (int a = 1; a < 6; a++) {
        N_y_a(0,a) = init_pop(a-1);
    }

    // estimate seine age comps, first year
    matrix<Type> seine_age_comp_est(nyr, nage);         // seine age composition
    vector<Type> seine_selected(nyr);
    // estimate fully selected seine catch to calculate age comp estimate
    seine_selected(0) = N_y_a.row(0).dot(seine_selex.matrix().transpose());      
    for (int a = 0; a < nage; a++) {
        seine_age_comp_est(0, a) = seine_selex(a)*N_y_a(0, a) / seine_selected(0);
    }


    // estimate seine catch,  first year
    vector<Type> seine_catch_est(nyr);
    seine_catch_est(0) = seine_yield(0) / (seine_age_comp_est.row(0).dot(waa.row(0))); 

    // total spring catch-at-age, first year
    matrix<Type> spring_removals(nyr, nage);
    spring_removals = spring_removals.setZero(nyr, nage);
    spring_removals.row(0) = seine_age_comp_est.row(0)*seine_catch_est(0) + pk*pound_catch.row(0) + gillnet_catch.row(0);

    // ---- model naa/caa dynamics 1981 to 1984 ---- //
    // Since the numbers-at-age matrix is only initialized up to age-5 (i.e. 
    // the log_initpop parameter is a vector of length 5), the older age classes 
    // have entries of 0. We will assume the largest initial age class is a plus
    // group. Hence, from 1981 to 1984 the numbers-at-age matrix will be 
    // calculated using incrementally larger plus groups until all age classes
    // are populated (i.e. age-5 is plus group in 1980, age-6 is plus group in 
    // 1981, age-7 is plus group in 1982, age-8 is plus group in 1983, age-9 
    // is plus group in 1984)

    int index_1984 = (1984-1980)+1;             // index of year where plus group becomes 9+
    Type naa_min = 0.01;
    Type naa_pen = 0.00;

    for (int y = 1; y < index_1984; y++) {
        int plus_group_index = 5 + y;                   // index of plus group 

        // pre-fishery numbers-at-age
        N_y_a(y,0) = age_0(y);    

        // calculate plus group catch (increments from age-6 to age-9+)
        Type spring_removals_sum = spring_removals.row(y-1).tail(nage-plus_group_index).sum();
        Type foodbait_sum = foodbait_catch.row(y-1).tail(nage-plus_group_index).sum();

        // populate age classes
        for (int a = 1; a < nage; a++) {
            if (a < plus_group_index) {
                // age-1 to plus group
                Type naa_pos = (((N_y_a(y-1,a-1) - spring_removals(y-1,a-1)) * summer_survival(y-1,a-1)) - foodbait_catch(y-1,a-1)) * winter_survival(y,a-1);
                if (naa_pos <= naa_min) {
                    N_y_a(y, a) = naa_min;
                    naa_pen += .01*pow(N_y_a(y, a) - naa_min, 2);
                } else if (naa_pos > naa_min) {
                    N_y_a(y, a) = naa_pos;
                }
            } else if (a == plus_group_index) {
                // plus group
                Type naa_pos = (((N_y_a(y-1,a-1) - spring_removals_sum) * summer_survival(y-1,a-1)) - foodbait_sum) * winter_survival(y,a-1);
                if (naa_pos <= naa_min) {
                    N_y_a(y, a) = naa_min;
                    naa_pen += .01*pow(naa_pos - naa_min, 2);
                } else if (naa_pos > naa_min) {
                    N_y_a(y, a) = naa_pos;
                }
            } else if (a > plus_group_index) {
                // numbers-at-age set to 0 if age greater than plus group
                N_y_a(y,a) = 0;
            }
        }

        // seine age comps estimates
        seine_selected(y) = N_y_a.row(y).dot(seine_selex.matrix().transpose());      
        for (int a = 0; a < nage; a++) {
            seine_age_comp_est(y, a) = seine_selex(a)*N_y_a(y, a) / seine_selected(y);
        }

        // estimate seine catch
        seine_catch_est(y) = seine_yield(y) / (seine_age_comp_est.row(y).dot(waa.row(y))); 

        // total spring catch-at-age
        spring_removals.row(y) = seine_age_comp_est.row(y)*seine_catch_est(y) + pk*pound_catch.row(y) + gillnet_catch.row(y);

    }

    // ---- model naa/caa dynamics 1985 to nyr ---- //
    // from 1985 on, the plus group is in a static age class (9+)

    for (int y = index_1984; y < nyr; y++) { 

        // pre-fishery numbers-at-age
        N_y_a(y,0) = age_0(y);    

        // populate age classes
        for (int a = 1; a < nage; a++) {
            Type naa_pos = (((N_y_a(y-1,a-1) - spring_removals(y-1,a-1)) * summer_survival(y-1,a-1)) - foodbait_catch(y-1,a-1)) * winter_survival(y,a-1);
            if (a == nage-1) {
                // plus group
                naa_pos += (((N_y_a(y-1,nage-1) - spring_removals(y-1,nage-1)) * summer_survival(y-1,nage-1)) - foodbait_catch(y-1,nage-1)) * winter_survival(y,nage-1);
            }
            if (naa_pos <= naa_min) {
                N_y_a(y, a) = naa_min;
                naa_pen += .01*pow(N_y_a(y, a) - naa_min, 2);
            } else if (naa_pos > naa_min) {
                N_y_a(y, a) = naa_pos;
            }
        }

        // seine age comps estimates
        seine_selected(y) = N_y_a.row(y).dot(seine_selex.matrix().transpose());      
        for (int a = 0; a < nage; a++) {
            seine_age_comp_est(y, a) = seine_selex(a)*N_y_a(y, a) / seine_selected(y);
        }

        // estimate seine catch
        seine_catch_est(y) = seine_yield(y) / (seine_age_comp_est.row(y).dot(waa.row(y))); 

        // total spring catch-at-age
        spring_removals.row(y) = seine_age_comp_est.row(y)*seine_catch_est(y) + pk*pound_catch.row(y) + gillnet_catch.row(y);

    }

    // ---- post-fishery spawning numbers-at-age matrix ---- //
    // calculate post-fishery spawning abundance 
    // spring removals are calculated here without incorporating assumed mortality rate
    // for pound spawn-on-kelp (pk), so the surviving herring that spawn in that fishery 
    // are not included in the post-fishery spawning abundance as their eggs are 
    // harvested 

    matrix<Type> Ntilde_y_a(nyr, nage);
    matrix<Type> spawn_removals(nyr, nage);
    spawn_removals = spawn_removals.setZero();
    Type ntilde_pen = 0.00;


    // ---- Penalize negative numbers-at-age ---- //

    // ---- 1981 to 1984 ---- //
    for (int y = 0; y < index_1984; y++) {
        int plus_group_index = 5 + y;                   // index of plus group

        // calculate plus group catch (increments from age-6 to age-9+)
        Type seine_sum = (seine_age_comp_est.row(y).tail(nage-plus_group_index)*seine_catch_est(y)).sum(); 
        Type pound_sum = pound_catch.row(y).tail(nage-plus_group_index).sum(); 
        Type gillnet_sum = gillnet_catch.row(y).tail(nage-plus_group_index).sum();

        for(int a = 0; a < nage; a++) {
            if (a < plus_group_index) {
                spawn_removals(y, a) = seine_age_comp_est(y,a)*seine_catch_est(y) + pound_catch(y,a) + gillnet_catch(y,a);
                Type ntilde_pos = maturity(a)*(N_y_a(y, a) - spawn_removals(y, a));
                if ((N_y_a(y, a) - spawn_removals(y, a)) <= naa_min) {
                    Ntilde_y_a(y, a) = maturity(a)*(spawn_removals(y, a) + naa_min);
                    ntilde_pen += .01*pow((N_y_a(y, a) - spawn_removals(y, a)) - naa_min, 2);
                } else if ((N_y_a(y, a) - spawn_removals(y, a)) > naa_min) {
                    Ntilde_y_a(y, a) = ntilde_pos;
                }
            } else if (a == plus_group_index) {
                spawn_removals(y, a) = seine_sum + pound_sum + gillnet_sum;
                Type ntilde_pos = maturity(a)*(N_y_a(y, a) - spawn_removals(y, a));
                if ((N_y_a(y, a) - spawn_removals(y, a)) <= naa_min) {
                    Ntilde_y_a(y, a) = maturity(a)*(spawn_removals(y, a) + naa_min);
                    ntilde_pen += .01*pow((N_y_a(y, a) - spawn_removals(y, a)) - naa_min, 2);                    
                } else if ((N_y_a(y, a) - spawn_removals(y, a)) > naa_min) {
                    Ntilde_y_a(y, a) = ntilde_pos;
                }
            } else if (a > plus_group_index) {
                spawn_removals(y, a) = 0;
                Ntilde_y_a(y, a) = 0;
            }
        }
    }

    // ---- 1985 to nyr ---- //
    for (int y = index_1984; y < nyr; y++) {
        for(int a = 0; a < nage; a++) {
            spawn_removals(y, a) = seine_age_comp_est(y,a)*seine_catch_est(y) + pound_catch(y,a) + gillnet_catch(y,a);
            Type ntilde_pos = maturity(a)*(N_y_a(y, a) - spawn_removals(y, a));
            if ((N_y_a(y, a) - spawn_removals(y, a)) <= naa_min) {
                Ntilde_y_a(y, a) = maturity(a)*(spawn_removals(y, a) + naa_min);
                ntilde_pen += .01*pow((N_y_a(y, a) - spawn_removals(y, a)) - naa_min, 2);
            } else if ((N_y_a(y, a) - spawn_removals(y, a)) > naa_min) {
                Ntilde_y_a(y, a) = ntilde_pos;
            }
        }
    }

    // ---- biomass estimates ---- //

    vector<Type> B_y(nyr);                              // pre-fishery total biomass
    matrix<Type> Btilde_y(nyr, nage);                   // pre-fishery spawning biomass
    vector<Type> Btilde_post_y(nyr);                    // post-fishery spawning biomass

    B_y = (N_y_a.array() * waa.array()).rowwise().sum();
    Btilde_y = (N_y_a.array() * waa.array()).matrix() * maturity.matrix();
    Btilde_post_y = (Ntilde_y_a.array() * waa.array()).rowwise().sum();

    // ---- biomass forecast ---- //

    vector<Type> projected_N_y_a(nage);                 // numbers-at-age forecast
    vector<Type> waa_5yr_mean(nage);                    // weight-at-age 5yr mean
    vector<Type> winter_survival_forecast(nage);        // winter survival for forecast
    Type Btilde_forecast;                               // spawning biomass forecast
    
    // set up winter survival forecast vector
    for (int a = 0; a < nage; a++) {
        if(a < nage - 1) {
            winter_survival_forecast(a) = exp(-0.5*Z_0_8);
        } else if (a == nage-1) { 
            winter_survival_forecast(a) = exp(-0.5*Z_9);            
        }
    }

    // project naa matrix forward 1 year
    for (int a = 0; a < nage; a++) {
        if (a < 4) {
            projected_N_y_a(a) = 0;   
        } else if (a >= 4) {
            projected_N_y_a(a) = ((N_y_a(nyr-1, a-1)-(seine_age_comp_est(nyr-1,a-1)*seine_catch_est(nyr-1)+gillnet_catch(nyr-1,a-1)+pk*pound_catch(nyr-1,a-1)))*summer_survival(nyr-1,a-1)-foodbait_catch(nyr-1,a-1))*winter_survival_forecast(a);
        } 
        if (a == nage-1) {
            // plus group
            projected_N_y_a(a) += ((N_y_a(nyr-1, a)-(seine_age_comp_est(nyr-1,a)*seine_catch_est(nyr-1)+gillnet_catch(nyr-1,a)+pk*pound_catch(nyr-1,a)))*summer_survival(nyr-1,a)-foodbait_catch(nyr-1,a))*winter_survival_forecast(a);
        }
    }


    // take 5-year mean of waa matrix
    for (int a = 0; a < nage; a++) {
        waa_5yr_mean(a) = waa.col(a).tail(5).sum() / 5;
    }

    // pre-fishery spawning biomass forecast
    for (int a = 0; a < nage; a++) {
        Btilde_forecast += maturity(a)*projected_N_y_a(a)*waa_5yr_mean(a);
    }

    // ----------- OTHER ESTIMATES IN LIKELIHOOD EXPRESSIONS ------------ //

    // ---- estimate spawning age compositions ---- //

    matrix<Type> spawn_age_comp_est(nyr, nage);         // spawning age composition
    vector<Type> total_spawners(nyr);                   // all spawners across ages
    total_spawners = N_y_a * maturity; 

    for (int y = 0; y < nyr; y++) {
        for (int a = 0; a < nage; a++) {
            spawn_age_comp_est(y, a) = maturity(a)*N_y_a(y, a) / total_spawners(y);
        }
    }

    // ---- estimate ADFG hydroacoustic biomass ---- //

    vector<Type> Hhat_adfg_y = Btilde_y * exp(adfg_hydro_q);

    // ---- estimate PWSSC hydroacoustic biomass ---- //

    vector<Type> Hhat_pwssc_y = Btilde_y * exp(pwssc_hydro_q);

    // ---- estimate naturally spawned eggs ---- //

    vector<Type> Ehat_y(nyr);
    vector<int> fecundity_years(nyr);
    fecundity_years = fecundity_years.setZero();

    // the following finds for which years the fecundity index exists 
    for (int y = 0; y < nyr; y++) {
        for (int a = 0; a < nage; a++) {
            fecundity_years(y) += (fecundity(y,a) != -9);
        }
        if (fecundity_years(y) > 0) fecundity_years(y) = 1; 
    }

    Ehat_y = Ehat_y.setZero();
    for (int y = 0; y < nyr; y++) {
        if (fecundity_years(y) == 0) continue;          // skip years with no fecundity index
        Ehat_y(y) = pow(10, -6) * perc_female(y) * Ntilde_y_a.row(y).dot(fecundity.row(y));
    }

    // ---- estimate mile-days milt ---- //

    vector<Type> That_y(nyr);
    for (int y = 0; y < nyr; y++) {
        That_y(y) = ((1 - perc_female(y)) * Btilde_post_y(y)) / exp(logmdm_c);
    }

    // ---- estimate juvenile aerial survey ---- //

    vector<Type> Jhat_y(nyr);

    for (int y = 0; y < nyr; y++){
        Jhat_y(y) = N_y_a(y,1) * exp(log_juvenile_q); 
    }

    // ------------------------- LIKELIHOODS ---------------------------- //

    // ---- L1: purse-seine age-composition ---- //

    // Calculate likelihoods from 1980 to 1984
    // increments plus group from age-5 in 1980 to age-9 in 1984 in a similar 
    // fashion to the numbers-at-age calculation 

    Type L1 = 0.0;
    vector<Type> L1_years(nyr);
    L1_years = L1_years.setZero();

    for (int y = 0; y < index_1984; y++) {
        if (seine_ess(y) == -9) continue;               // skip years with no fishery
        int plus_group_index = 5 + y;                   // index of plus group 
        Type seine_plus_group = seine_age_comp.row(y).tail(nage-plus_group_index).sum();
        for (int a = 3; a < plus_group_index; a++) {
            // protects observed age classes with 0% from blowing up likelihood
            if (seine_age_comp(y,a) <= 0) continue; 
            L1_years(y) += seine_age_comp(y,a)*log(seine_age_comp_est(y,a)/seine_age_comp(y,a));
        }
        L1_years(y) += seine_plus_group*log(seine_age_comp_est(y,plus_group_index)/seine_plus_group);
        L1 -= seine_ess(y)*L1_years(y);
    } 

    // Calculate likelihoods from 1984 to nyr
    
    for (int y = index_1984; y < nyr; y++) {
        if (seine_ess(y) == -9) continue;              // skip years with no fishery
        for (int a = 3; a < nage; a++) {
            if (seine_age_comp(y,a) <= 0) continue;
            L1_years(y) += seine_age_comp(y,a)*log(seine_age_comp_est(y,a)/seine_age_comp(y,a));
        }
        L1 -= seine_ess(y)*L1_years(y);
    } 

    // ---- L2: spawn survey age-composition ---- //

    // Calculate likelihoods from 1980 to 1984
    // increments plus group from age-5 in 1980 to age-9 in 1984 in a similar 
    // fashion to the numbers-at-age calculation 

    Type L2 = 0.0;
    vector<Type> L2_years(nyr);
    L2_years = L2_years.setZero();

    for (int y = 0; y < index_1984; y++) {
        if (spawn_ess(y) == -9) continue;               // skip years with no survey
        int plus_group_index = 5 + y;                   // index of plus group 
        Type spawn_plus_group = spawn_age_comp.row(y).tail(nage-plus_group_index).sum();
        for (int a = 3; a < plus_group_index; a++) {
            // protects observed age classes with 0% from blowing up likelihood
            if (spawn_age_comp(y,a) <= 0) continue;
            L2_years(y) += spawn_age_comp(y,a)*log(spawn_age_comp_est(y,a)/spawn_age_comp(y,a));
        }
        L2_years(y) += spawn_plus_group*log(spawn_age_comp_est(y,plus_group_index)/spawn_plus_group);
        L2 -= spawn_ess(y)*L2_years(y);
    } 

    // Calculate likelihoods from 1984 to nyr
    
    for (int y = index_1984; y < nyr; y++) {
        if (spawn_ess(y) == -9) continue;               // skip years with no survey
        for (int a = 3; a < nage; a++) {
            if (spawn_age_comp(y,a) <= 0) continue;
            L2_years(y) += spawn_age_comp(y,a)*log(spawn_age_comp_est(y,a)/spawn_age_comp(y,a));
        }
        L2 -= spawn_ess(y)*L2_years(y);
    } 

    // ---- L3: Number of eggs deposited ---- //

    Type L3 = 0.0;
    vector<Type> L3_total_var(nyr);
    L3_total_var = L3_total_var.setZero();


    // number of years in egg index
    int n_egg = 0;

    L3_total_var = L3_total_var.setZero();
    for (int y = 0; y < nyr; y++) {
        if (egg(y) == -9) continue;                     // skip years with no egg index
        L3_total_var(y) = pow(egg_se(y), 2) + pow(egg_add, 2);
        L3 += log(pow(L3_total_var(y), .5));
        L3 += (pow(log(Ehat_y(y)/egg(y)), 2) / (2*L3_total_var(y)));
        n_egg += 1; 
    }
    // L3 *= n_egg;

    // ---- L4: ADFG Hydroacoustic biomass ---- //

    Type L4 = 0.0;

    // number of years in ADFG hydroacoustic survey
    int n_adfg_hydro = 0;

    for (int y = 0; y < nyr; y++) {
        if (adfg_hydro(y) == -9) continue;              // skip years when no ADFG hydro survey
        n_adfg_hydro += 1;
        L4 += pow(log(Hhat_adfg_y(y)/adfg_hydro(y)), 2);
    }
    L4 /= 2*pow(adfg_hydro_add_var, 2);
    L4 += n_adfg_hydro*log(adfg_hydro_add_var);

    // ---- L5: PWSSC Hydroacoustic biomass ---- //

    Type L5 = 0.0;
    vector<Type> L5_total_var(nyr);
    L5_total_var = L5_total_var.setZero();

    // number of years in PWSSC hydroacoustic survey
    int n_pwssc_hydro = 0;

    for (int y = 0; y < nyr; y++) {
        if (pwssc_hydro_se(y) == -9) continue;          // skip years when no PWSSC hydro survey
        L5_total_var(y) = pow(pwssc_hydro_se(y), 2) + pow(pwssc_hydro_add_var, 2);
        L5 += log(pow(L5_total_var(y), .5));
        L5 += (pow(log(Hhat_pwssc_y(y)/pwssc_hydro(y)), 2) / (2*L5_total_var(y)));
        n_pwssc_hydro += 1;
    }
    // L5 *= n_pwssc_hydro;

    // ---- L6: Mile-days milt ---- //

    Type L6 = 0.0;

    // number of years in MDM index
    int n_mdm = 0;
    
    for (int y = 0; y < nyr; y++) {
        if (mdm(y) == -9) continue;
        n_mdm += 1;
        L6 += pow(log(That_y(y) / mdm(y)), 2);
    }
    L6 /= 2*pow(milt_add_var, 2);
    L6 += n_mdm*log(milt_add_var);

    // ---- L7: Juvenile aerial survey ---- //

    Type L7 = 0.0;
    vector<Type> L7_var(nyr);

    for (int y = 0; y < nyr; y++) {
        L7_var(y) = Jhat_y(y) + (pow(Jhat_y(y), 2)/juvenile_overdispersion);
        if (juvenile_survey(y) == -9) continue;
        L7 -= dnbinom2(Type(juvenile_survey(y)), Jhat_y(y), L7_var(y), true);
    }

    // ---- total likelihood ---- //

    Type negLogLik = 0.0;                               
    negLogLik = L1 + L2 + L3 + L4 + L5 + L6 + L7 + 1000*naa_pen + 1000*ntilde_pen + 1000*winter_surv_pen + 1000*summer_surv_pen;

    // ---- incorporate prior densities ---- //

    for (int y = 0; y < nyr-1; y++) {
        negLogLik -= dunif(annual_age0devs(y), Type(-10), Type(10), true);
    }
    for (int y = 0; y < 5; y++) {
        negLogLik -= dunif(loginit_pop(y), Type(3.00), Type(8.00), true);
    }
    // negLogLik -= dnorm(MeanAge0, Type(0), Type(sigma_age0devs), true);
    // negLogLik -= dunif(sigma_age0devs, Type(0.00), Type(2.0), true);
    negLogLik -= dunif(log_juvenile_q, Type(-5), Type(8.0), true);
    negLogLik -= dunif(Z_9, Type(0.3), Type(1.6), true);
    for (int b = 0; b < n_covs; b++) {
        negLogLik -= dunif(beta_mortality(b), Type(-10.00), Type(10.00), true);
    }
    negLogLik -= dunif(logmdm_c, Type(2.3), Type(9.0), true);
    negLogLik -= dunif(adfg_hydro_q , Type(-5), Type(5), true);
    negLogLik -= dunif(pwssc_hydro_q, Type(-5), Type(5), true);
    negLogLik -= dunif(VHSV_age3_4_mort_93, Type(0.0), Type(5.0), true);
    negLogLik -= dunif(ICH_age5_8_mort_93, Type(0.0), Type(5.0), true);
    negLogLik -= dbeta(mat_age3, Type(9.0), Type(11.0), true);
    negLogLik -= dbeta(mat_age4, Type(18.0), Type(2.0), true);
    negLogLik -= dunif(seine_selex_alpha, Type(3.0), Type(5.0), true);
    negLogLik -= dunif(seine_selex_beta, Type(1.0), Type(7.0), true);
    negLogLik -= dnorm(milt_add_var, Type(0.33), Type(0.10), true);
    negLogLik -= dnorm(adfg_hydro_add_var, Type(0.30), Type(0.08), true);
    negLogLik -= dnorm(pwssc_hydro_add_var, Type(0.32), Type(0.08), true);
    negLogLik -= dunif(juvenile_overdispersion, Type(0.01), Type(4.00), true);


    // ------------------------------------------------------------------ //
    //                         Report Section                             //
    // ------------------------------------------------------------------ //

    REPORT(seine_age_comp_est);
    REPORT(spawn_age_comp_est);
    REPORT(age_0);
    REPORT(N_y_a);
    REPORT(Ntilde_y_a);
    REPORT(B_y);
    REPORT(Btilde_y);
    REPORT(Btilde_post_y);
    REPORT(spawn_removals);
    REPORT(summer_survival);
    REPORT(winter_survival);
    REPORT(L1);
    REPORT(L2);
    REPORT(L3);
    REPORT(L4);
    REPORT(L5);
    REPORT(L6);
    REPORT(L7);
    REPORT(negLogLik);
    REPORT(naa_pen);
    REPORT(ntilde_pen);
    REPORT(Ehat_y);
    REPORT(L3_total_var);
    REPORT(Hhat_adfg_y);
    REPORT(Hhat_pwssc_y);
    REPORT(L5_total_var);
    REPORT(That_y);
    REPORT(Jhat_y);
    REPORT(L7_var);
    REPORT(projected_N_y_a);
    REPORT(waa_5yr_mean);
    REPORT(Btilde_forecast);    
    
    return negLogLik;
  
}
