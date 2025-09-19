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
// C++ program file:  PWS_ASA.cpp                                             //
// R control file:  run_basa.r                                                //
//                                                                            //
// Input Data files:                                                          //
//        PWS_ASA.dat:             Model input (surveys, catches, etc.)       //
//                                                                            //
// Output files:                                                              //
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

template <class Type>
Type posfun(Type x, Type eps, Type &pen){
    pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
    return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template <class Type>
Type calculate_logistic(int a, Type a50, Type a95){
    Type logistic = 0.0;
    logistic = 1.0 / (1.0 + exp(-log(19)*(a - a50)/(a95 - a50)));
    return logistic;
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
    DATA_IVECTOR(adfg_hydro);               // ADFG hydroacoustic survey index    
    DATA_IVECTOR(pwssc_hydro);              // PWSSC hydroacoustic survey index
    DATA_VECTOR(pwssc_hydro_se);            // PWSSC hydroacoustic survey se's
    DATA_MATRIX(seine_age_comp);             // purse-seine age compositions
    DATA_MATRIX(spawn_age_comp);             // spawn survey age compositions
    DATA_IVECTOR(juvenile_survey);          // aerial juvenile school survey index

    // ---- model/PWS_ASA(covariate).ctl ---- //
    DATA_MATRIX(disease_covs);              // excess mortality due to disease
    DATA_IMATRIX(mort_age_impact);           // which age does the mort cov impact?

    // ---- model/PWS_ASA_disease.dat ---- //
    DATA_MATRIX(vhsv_obs);       // VHSV seroprevalence
    DATA_INTEGER(vhsv_est_start);            // start index seroprevalence model estimation
    DATA_INTEGER(vhsv_obs_start);            // start index seroprevalence time series

    // ---- model/agecomp_samp_sizes.txt ---- //
    DATA_IVECTOR(seine_sample_size);        // purse seine raw sample sizes
    DATA_IVECTOR(spawn_sample_size);        // spawn survey raw sample sizes
    DATA_IVECTOR(vhsv_sample_size);         // antibody raw sample sizes

    // ---- forecast controls ---- //
    DATA_INTEGER(recruitment_average_years);    // # years to average recruitments in forecast
    DATA_INTEGER(waa_average_years);            // # years to average WAA in forecast 
    DATA_INTEGER(disease_cov_average_years);    // # years to average disease covariates in forecast
    DATA_SCALAR(expected_spring_harvest);       // metric tons of expected harvest for milt forecast
    DATA_INTEGER(perc_female_forecast_years);   // # years to average percent females in milt forecast

    // ------------------------------------------------------------------ //
    //                        Parameter Section                           //
    // ------------------------------------------------------------------ //

    // ---- fixed parameters ---- //
    PARAMETER(pk);                      // pound spawn-on-kelp mortality proportion
    PARAMETER(egg_add);                 // egg deposition additional variance term
    PARAMETER(Z_0_8);                   // background mortality, ages 0-8
    // PARAMETER(sigma_age0devs);          // variance term for recruitment deviates
    
    // ---- estimated parameters ---- //
    PARAMETER(log_MeanAge0);            // mean recruits, log-space
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
    
    PARAMETER(vhs_inf_a50);
    PARAMETER(vhs_inf_delta);
    // PARAMETER(vhs_samp_a50);
    // PARAMETER(vhs_samp_a95);
    PARAMETER_VECTOR(vhs_inf_prob);
    PARAMETER_VECTOR(vhs_rec_prob);
    
    // ------------------------------------------------------------------ //
    //                        Procedure Section                           //
    // ------------------------------------------------------------------ //

    // ---- convert some parameters from log-space ---- //
    vector<Type> init_pop = exp(loginit_pop);           // initial population vector

    // ---- penalty objects ---- //
    int penCount = 0;                                   // count instances of penalties

    // ---- dummy objects ---- //
    Type dummy = 0;                                     // scalar for saving dummy variables for debugging
    vector<Type> dummy_vector(nyr);                     // dummy vector
    dummy_vector.setZero();
    matrix<Type> dummy_matrix(nyr, nage);               // dummy matrix
    dummy_matrix.setZero();

    // --------------------- POPULATION DYNAMICS ------------------------ //

    // ---- survival from natural mortality ---- //
    // calculate half-year survival using background mortality, plus-group
    // background mortality, excess mortality in 1992-1993 for the different
    // age groups, and 3 disease covariates

    // 2 matrices (summer and winter) will be used to store half-year survival 
    matrix<Type> summer_survival(nyr, nage);
    summer_survival.setZero();
    matrix<Type> summer_mortality_effect(nyr, nage);
    summer_mortality_effect.setZero(); 
    
    matrix<Type> winter_survival(nyr, nage);
    winter_survival.setZero();
    matrix<Type> winter_mortality_effect(nyr, nage);
    winter_mortality_effect.setZero(); 
    
    int index_1992 = 1992-1980;                 // index for 1992
    int n_covs = beta_mortality.size();
    
    matrix<Type> disease_covs_calc(nyr, n_covs);
    disease_covs_calc.setZero();

    Type min_mortality = .01;
    Type winter_surv_pen = 0.0;
    Type summer_surv_pen = 0.0;

    // calculate disease mortality effects
    for (int y = 0; y < nyr; y++) {
        for(int a = 0; a < nage; a++) {

            for(int b = 0; b < n_covs; b++) {

                if (disease_covs(y, b) == -9) {
                    disease_covs_calc(y, b) = 0;
                } else if (disease_covs(y, b) != -9) {
                    disease_covs_calc(y, b) = disease_covs(y, b);
                }

                // temp measure to fix vhsv beta mortality param to zero
                if (b == 2) {
                    disease_covs_calc(y, b) = 0;
                }

                summer_mortality_effect(y, a) += beta_mortality(b)*disease_covs_calc(y, b)*mort_age_impact(a, b);
                if(y < nyr-1) {
                    winter_mortality_effect(y+1, a) += beta_mortality(b)*disease_covs_calc(y, b)*mort_age_impact(a, b);
                }              
            }

        }
    }

    // calculate survival matrices
    for (int y = 0; y < nyr; y++) {
        for (int a = 0; a < (nage-1); a++) {
            // ages 0-8
            summer_survival(y, a) = exp(-0.5*Z_0_8 - summer_mortality_effect(y, a));
            winter_survival(y, a) = exp(-0.5*Z_0_8 - winter_mortality_effect(y, a));

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
            
            // penalize survivals >= 1
            summer_survival(y, a) = 1 - posfun(1-summer_survival(y, a), min_mortality, summer_surv_pen);
            if(summer_survival(y, a) >= 1 - min_mortality) penCount += 1;

            winter_survival(y, a) = 1 - posfun(1-winter_survival(y, a), min_mortality, winter_surv_pen);
            if(winter_survival(y, a) >= 1 - min_mortality) penCount += 1;

        } 

        // plus group (9+)
        if (y == 0) {

            summer_survival(y, nage-1) = exp(-0.5*Z_9 - summer_mortality_effect(y, nage-1));
            winter_survival(y, nage-1) = exp(-0.5*Z_9 - winter_mortality_effect(y, nage-1));
            
        } else if (y > 0) {
            
            summer_survival(y, nage-1) = summer_survival(y-1, nage-1) * (summer_survival(y, nage-2)/summer_survival(y-1, nage-2));
            winter_survival(y, nage-1) = winter_survival(y-1, nage-1) * (winter_survival(y, nage-2)/winter_survival(y-1, nage-2));

        }
        
        // penalize survivals >= 1 for plus group
        summer_survival(y, nage-1) = 1 - posfun(1-summer_survival(y, nage-1), min_mortality, summer_surv_pen);
        if(summer_survival(y, nage-1) >= 1 - min_mortality) penCount += 1;
        
        winter_survival(y, nage-1) = 1 - posfun(1-winter_survival(y, nage-1), min_mortality, winter_surv_pen);
        if(winter_survival(y, nage-1) >= 1 - min_mortality) penCount += 1;
        
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



    // ---- VHSV seroprevalence dynamics ---- //
    
    // calculate VHS infection vulnerability
    // the proportion of each age class vulnerable to transmissions (i.e. well mixed with the infected population)  
    // matrix<Type> disease_inf_vul(nyr, nage);
    // disease_inf_vul.setZero();
    // for (int y = 0; y < nyr; y++){
    //     for (int a = 0; a < nage; a++){
    //         disease_inf_vul(y, a) = calculate_logistic(a, vhs_inf_a50, vhs_inf_a95);
    //     }
    // }
    Type vhs_inf_a95 = vhs_inf_a50 + vhs_inf_delta;
    vector<Type> disease_inf_vul(nage);
    disease_inf_vul.setZero();
    for (int a = 0; a < nage; a++){
        disease_inf_vul(a) = calculate_logistic(a, vhs_inf_a50, vhs_inf_a95);
    }

    
    // calculate VHS survey selectivity
    // the proportion of infected fish each age class vulnerable to detection in survey 
    // assumes survey sample vulnerability = maturity
    // vector<Type> disease_samp_vul(nage);
    // disease_samp_vul.setZero();
    // disease_samp_vul = maturity;
    // for (int y = 0; y < nyr; y++){
    //     for (int a = 0; a < nage; a++){
    //         disease_samp_vul(y, a) = calculate_logistic(a, vhs_samp_a50, vhs_samp_a95);
    //     }
    // }

    vector<Type> disease_samp_vul(nage);
    disease_samp_vul.setZero();
    disease_samp_vul = maturity;
    // for (int a = 0; a < nage; a++){
    //     disease_samp_vul(a) = calculate_logistic(a, vhs_samp_a50, vhs_samp_a95);
    // }

    
    // these are vectors constructed from the estimated infection/recovery 
    // parameters with the same length as the model time series for indexing purposes
    
    vector<Type> vhs_inf_prob_y(nyr);
    vhs_inf_prob_y.setZero();
    vector<Type> vhs_rec_prob_y(nyr);
    vhs_rec_prob_y.setZero();
    
    for (int y = 0; y < nyr; y++) {
        
        int vhs_index = y - vhsv_est_start;
        if (vhs_index >= 0) {
            vhs_inf_prob_y(y-1) = vhs_inf_prob(vhs_index);
            vhs_rec_prob_y(y-1) = vhs_rec_prob(vhs_index);
        }

    }
    
    // ---- recruitment ---- //
    vector<Type> age_0(nyr);
    age_0.setZero();

    for (int y = 0; y < nyr-1; y++) {
        age_0(y) = exp(log_MeanAge0 + annual_age0devs(y));
    }

    // set the recruitment deviation to 0 in the final year  
    age_0(nyr-1) = exp(log_MeanAge0);                                       


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
    seine_age_comp_est.setZero();
    vector<Type> seine_selected(nyr);
    seine_selected.setZero();
    
    // estimate fully selected seine catch to calculate age comp estimate
    seine_selected(0) = N_y_a.row(0).dot(seine_selex.matrix().transpose());      
    for (int a = 0; a < nage; a++) {
        seine_age_comp_est(0, a) = seine_selex(a)*N_y_a(0, a) / seine_selected(0);
    }

    // estimate seine catch,  first year
    vector<Type> seine_catch_est(nyr);
    seine_catch_est.setZero();

    seine_catch_est(0) = seine_yield(0) / (seine_age_comp_est.row(0).dot(waa.row(0))); 

    // total spring catch-at-age, first year
    matrix<Type> spring_removals(nyr, nage);
    spring_removals.setZero(nyr, nage);

    spring_removals.row(0) += seine_age_comp_est.row(0)*seine_catch_est(0); 
    spring_removals.row(0) += pk*pound_catch.row(0); 
    spring_removals.row(0) += gillnet_catch.row(0);

    // ---- model naa/caa dynamics 1981 to 1984 ---- //
    // Since the numbers-at-age matrix is only initialized up to age-5 (i.e. 
    // the log_initpop parameter is a vector of length 5), the older age classes 
    // have entries of 0. We will assume the largest initial age class is a plus
    // group. Hence, from 1981 to 1984 the numbers-at-age matrix will be 
    // calculated using incrementally larger plus groups until all age classes
    // are populated (i.e. age-5 is plus group in 1980, age-6 is plus group in 
    // 1981, age-7 is plus group in 1982, age-8 is plus group in 1983, age-9 
    // is plus group in 1984)

    // vhs survival ignored since there is no seroprevalence data for these years

    int index_1984 = (1984-1980)+1;             // index of year where plus group becomes 9+
    Type naa_min = 0.01;
    Type naa_pen = 0.00;

    for (int y = 1; y < index_1984; y++) {
        int plus_group_index = 5 + y;                   // index of plus group 

        // pre-fishery numbers-at-age
        N_y_a(y,0) = age_0(y);    

        // calculate plus group catch (increments from age-6 to age-9+)
        Type spring_removals_sum = spring_removals.row(y-1).tail(nage-plus_group_index+1).sum();
        Type foodbait_sum = foodbait_catch.row(y-1).tail(nage-plus_group_index+1).sum();

        // populate age classes
        for (int a = 1; a < nage; a++) {
            if (a < plus_group_index) {                
                
                // age-1 to plus group
                N_y_a(y, a) = N_y_a(y-1, a-1) - spring_removals(y-1, a-1);
                N_y_a(y, a) *= summer_survival(y-1, a-1);
                N_y_a(y, a) -= foodbait_catch(y-1, a-1);
                N_y_a(y, a) *= winter_survival(y, a-1);
                      
                // penalize negative naa
                N_y_a(y, a) = posfun(N_y_a(y, a), naa_min, naa_pen);
                if(N_y_a(y, a) <= naa_min) penCount += 1;
                

            } else if (a == plus_group_index) {

                // plus group
                N_y_a(y, a) = N_y_a(y-1, a-1) - spring_removals_sum;
                N_y_a(y, a) *= summer_survival(y-1, a-1);
                N_y_a(y, a) -= foodbait_sum;
                N_y_a(y, a) *= winter_survival(y, a-1);

                // penalize negative naa
                N_y_a(y, a) = posfun(N_y_a(y, a), naa_min, naa_pen);
                if(N_y_a(y, a) <= naa_min) penCount += 1;

            } else if (a > plus_group_index) {

                // numbers-at-age set to 0 if age greater than plus group
                N_y_a(y, a) = 0;
            
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
        spring_removals.row(y) = seine_age_comp_est.row(y)*seine_catch_est(y); 
        spring_removals.row(y) += pk*pound_catch.row(y); 
        spring_removals.row(y) += gillnet_catch.row(y);

    }

    // ---- model naa/caa dynamics 1985 to nyr ---- //
    // from 1985 on, the plus group is in a static age class (9+)
    // VHS seroprevalence model dynamics begin

    // objects for seroprevalence model
    vector<Type> N_pre_vhs(nage);
    N_pre_vhs = N_pre_vhs.setZero(nage);

    matrix<Type> vhs_survival(nyr, nage);
    vhs_survival.setConstant(1);
    
    matrix<Type> vhs_immune(nyr, nage);
    vhs_immune.setZero();
    
    matrix<Type> vhs_susceptible(nyr, nage);
    vhs_susceptible.setConstant(1);

    for (int y = index_1984; y < nyr; y++) { 

        // pre-fishery numbers-at-age
        N_y_a(y,0) = age_0(y);    
        N_pre_vhs(0) = age_0(y);    

        // populate older age classes
        for (int a = 1; a < nage; a++) {
            
            // vhs survival ages 0 through 8
            vhs_survival(y-1, a-1) = 1 - disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1);
            vhs_survival(y-1, a-1) += disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1);
            
            // vhs immunity/susceptibility ages 1 through 8
            if (a < nage-1) {
                
                vhs_immune(y, a) = vhs_immune(y-1, a-1); 
                vhs_immune(y, a) += disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1);
                vhs_immune(y, a) /= (1 - disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1)) + (disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1));
                
                vhs_susceptible(y, a) = 1 - vhs_immune(y, a);
                
            }
            
            // pre-VHS mortality numbers-at-age
            N_pre_vhs(a) = N_y_a(y-1, a-1) - spring_removals(y-1, a-1);
            N_pre_vhs(a) *= summer_survival(y-1, a-1);
            N_pre_vhs(a) -= foodbait_catch(y-1, a-1);
            N_pre_vhs(a) *= winter_survival(y, a-1);
            
            // pre-fishery mortality numbers-at-age for start of year
            N_y_a(y, a) = N_pre_vhs(a)*vhs_survival(y-1, a-1);
            
            if (a == nage-1) {
                
                // VHS survival plus group
                vhs_survival(y-1, a) = 1 - vhs_susceptible(y-1, a)*vhs_inf_prob_y(y-1);
                vhs_survival(y-1, a) += vhs_susceptible(y-1, a)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1);
                
                // temp scalars for calculating VHS immunity/susceptibility plus groups
                // immune proportion of age-9's
                Type vhs_immune_age9 = 0.0;
                vhs_immune_age9 = vhs_immune(y-1, a-1);
                vhs_immune_age9 += disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1);
                vhs_immune_age9 /= (1 - disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1)) + (disease_inf_vul(a-1)*vhs_susceptible(y-1, a-1)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1)); 
                
                // immune proportion of age-10+
                Type vhs_immune_age10_plus = 0.0;
                vhs_immune_age10_plus = vhs_immune(y-1, a);
                vhs_immune_age10_plus += disease_inf_vul(a)*vhs_susceptible(y-1, a)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1);
                vhs_immune_age10_plus /= (1 - disease_inf_vul(a)*vhs_susceptible(y-1, a)*vhs_inf_prob_y(y-1)) + (disease_inf_vul(a)*vhs_susceptible(y-1, a)*vhs_inf_prob_y(y-1)*vhs_rec_prob_y(y-1));
                
                // numbers-at-age age-10+
                Type N_y_a_plus_group = 0.0;
                N_y_a_plus_group = N_y_a(y-1, a) - spring_removals(y-1, a);
                N_y_a_plus_group *= summer_survival(y-1, a);
                N_y_a_plus_group -= foodbait_catch(y-1, a);
                N_y_a_plus_group *= winter_survival(y, a);

                // fill in VHS immunity/susceptibility plus groups
                vhs_immune(y, a) = vhs_immune_age9*N_pre_vhs(a); 
                vhs_immune(y, a) += vhs_immune_age10_plus*N_y_a_plus_group;
                vhs_immune(y, a) /= N_pre_vhs(a) + N_y_a_plus_group;
                vhs_susceptible(y, a) = 1 - vhs_immune(y, a);
                N_pre_vhs(a) += N_y_a_plus_group;
                
                // accumulate age-10+ fish that survive VHS in numbers-at-age plus group 
                N_y_a(y, a) += N_y_a_plus_group*vhs_survival(y-1, a);
                
            }
            
            // penalize negative naa
            N_y_a(y, a) = posfun(N_y_a(y, a), naa_min, naa_pen);
            if(N_y_a(y, a) <= naa_min) penCount += 1;
            
        }

        // seine age comps estimates
        seine_selected(y) = N_y_a.row(y).dot(seine_selex.matrix().transpose());      
        for (int a = 0; a < nage; a++) {
            seine_age_comp_est(y, a) = seine_selex(a)*N_y_a(y, a) / seine_selected(y);
        }

        // estimate seine catch
        seine_catch_est(y) = seine_yield(y) / (seine_age_comp_est.row(y).dot(waa.row(y))); 

        // total spring catch-at-age
        spring_removals.row(y) = seine_age_comp_est.row(y)*seine_catch_est(y); 
        spring_removals.row(y) += pk*pound_catch.row(y); 
        spring_removals.row(y) += gillnet_catch.row(y);

    }

    // ---- post-fishery spawning numbers-at-age matrix ---- //
    // calculate post-fishery spawning abundance 
    // spring removals are calculated here without incorporating assumed mortality rate
    // for pound spawn-on-kelp (pk), so the surviving herring that spawn in that fishery 
    // are not included in the post-fishery spawning abundance as their eggs are 
    // harvested 

    matrix<Type> Ntilde_y_a(nyr, nage);
    Ntilde_y_a.setZero();

    matrix<Type> spawn_removals(nyr, nage);
    spawn_removals.setZero();

    Type ntilde_pen = 0.00;

    // ---- 1980 to 1984 ---- //
    for (int y = 0; y < index_1984; y++) {
        int plus_group_index = 5 + y;                   // index of plus group

        // calculate plus group catch (increments from age-6 to age-9+)
        Type seine_sum = (seine_age_comp_est.row(y).tail(nage-plus_group_index)*seine_catch_est(y)).sum(); 
        Type pound_sum = pound_catch.row(y).tail(nage-plus_group_index).sum(); 
        Type gillnet_sum = gillnet_catch.row(y).tail(nage-plus_group_index).sum();

        for(int a = 3; a < nage; a++) {
            if (a < plus_group_index) {

                spawn_removals(y, a) = seine_age_comp_est(y, a)*seine_catch_est(y);
                spawn_removals(y, a) += pound_catch(y, a); 
                spawn_removals(y, a) += gillnet_catch(y, a);
                
                Ntilde_y_a(y, a) = maturity(a)*(N_y_a(y, a) - spawn_removals(y, a));

                // penalize negative naa
                Ntilde_y_a(y, a) = posfun(Ntilde_y_a(y, a), naa_min, ntilde_pen);
                if(Ntilde_y_a(y, a) <= naa_min) penCount += 1;

            } else if (a == plus_group_index) {

                spawn_removals(y, a) = seine_sum + pound_sum + gillnet_sum;
                Ntilde_y_a(y, a) = maturity(a)*(N_y_a(y, a) - spawn_removals(y, a));

                // penalize negative naa
                Ntilde_y_a(y, a) = posfun(Ntilde_y_a(y, a), naa_min, ntilde_pen);
                if(Ntilde_y_a(y, a) <= naa_min) penCount += 1;

            } else if (a > plus_group_index) {

                spawn_removals(y, a) = 0;
                Ntilde_y_a(y, a) = 0;

            }
        }
    }

    // ---- 1985 to nyr ---- //
    for (int y = index_1984; y < nyr; y++) {
        for(int a = 3; a < nage; a++) {

            spawn_removals(y, a) = seine_age_comp_est(y, a)*seine_catch_est(y);
            spawn_removals(y, a) += pound_catch(y, a); 
            spawn_removals(y, a) += gillnet_catch(y, a);
            Ntilde_y_a(y, a) = maturity(a)*(N_y_a(y, a) - spawn_removals(y, a));

            // penalize negative naa
            Ntilde_y_a(y, a) = posfun(Ntilde_y_a(y, a), naa_min, ntilde_pen);
            if(Ntilde_y_a(y, a) <= naa_min) penCount += 1;

        }
    }


    // ---- biomass estimates ---- //
    
    vector<Type> B_y(nyr);                              // pre-fishery total biomass
    B_y.setZero();

    matrix<Type> Btilde_y(nyr, nage);                   // pre-fishery spawning biomass
    Btilde_y.setZero();

    vector<Type> Btilde_post_y(nyr);                    // post-fishery spawning biomass
    Btilde_post_y.setZero();

    B_y = (N_y_a.array() * waa.array()).rowwise().sum();
    Btilde_y = (N_y_a.array() * waa.array()).matrix() * maturity.matrix();
    Btilde_post_y = (Ntilde_y_a.array() * waa.array()).rowwise().sum();


    // derived vectors for plotting seroprevalence fits

    vector<Type> incidence_sp(nyr);
    incidence_sp.setZero();
    
    vector<Type> fatalities_sp(nyr);
    fatalities_sp.setZero();
    
    vector<Type> seroprev_sp(nyr);
    seroprev_sp.setZero();

    for(int y = 0; y < nyr; y++){

        if (y >= vhsv_est_start) {

            vector<Type> incidence_sp_numbers(nage);
            vector<Type> fatalities_sp_numbers(nage);
            vector<Type> seroprev_sp_numbers(nage);

            for(int a = 0; a < nage; a++){
                incidence_sp_numbers(a) = maturity(a)*N_y_a(y, a)*disease_inf_vul(a)*vhs_susceptible(y, a)*vhs_inf_prob_y(y);
                fatalities_sp_numbers(a) = maturity(a)*N_y_a(y, a)*disease_inf_vul(a)*vhs_susceptible(y, a)*vhs_inf_prob_y(y)*(1-vhs_rec_prob_y(y));
                seroprev_sp_numbers(a) = maturity(a)*N_y_a(y, a)*vhs_immune(y, a);
            }
                    
            incidence_sp(y) = incidence_sp_numbers.sum() / Ntilde_y_a.row(y).sum();
            fatalities_sp(y) = fatalities_sp_numbers.sum() / Ntilde_y_a.row(y).sum();
            seroprev_sp(y) = seroprev_sp_numbers.sum() / Ntilde_y_a.row(y).sum();
            
        }
    }

    // ----------------- ESTIMATED SURVEY QUANTITIES -------------------- //
    
    // ---- estimate spawning age compositions ---- //
    
    matrix<Type> spawn_age_comp_est(nyr, nage);         // spawning age composition
    spawn_age_comp_est.setZero();
    vector<Type> total_spawners(nyr);                   // all spawners across ages
    total_spawners.setZero();
    
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
    Ehat_y.setZero();

    vector<int> fecundity_years(nyr);
    fecundity_years.setZero();

    // the following finds for which years the fecundity index exists 
    for (int y = 0; y < nyr; y++) {
        for (int a = 0; a < nage; a++) {
            fecundity_years(y) += (fecundity(y, a) != -9);
        }
        if (fecundity_years(y) > 0) fecundity_years(y) = 1; 
    }

    for (int y = 0; y < nyr; y++) {
        if (fecundity_years(y) == 0) continue;          // skip years with no fecundity index
        Ehat_y(y) = pow(10, -6) * perc_female(y) * Ntilde_y_a.row(y).dot(fecundity.row(y));
    }

    // ---- estimate mile-days milt ---- //

    vector<Type> That_y(nyr);
    That_y.setZero();

    for (int y = 0; y < nyr; y++) {
        That_y(y) = ((1 - perc_female(y)) * Btilde_post_y(y)) / exp(logmdm_c);
    }

    // ---- estimate juvenile aerial survey ---- //
    
    vector<Type> Jhat_y(nyr);
    Jhat_y.setZero();
    
    for (int y = 0; y < nyr; y++){
        Jhat_y(y) = N_y_a(y,1) * exp(log_juvenile_q); 
    }

    // ---- estimate seroprevalence survey ---- //

    matrix<Type> vhsv_pred_numbers(nyr, nage*2);
    vhsv_pred_numbers.setZero();

    matrix<Type> vhsv_pred(nyr, nage*2);
    vhsv_pred.setZero();

    // the two age indices are used to alternate positive/negative seroprevalence
    for (int y = vhsv_est_start; y < nyr; y++) {

        for (int a = 0; a < nage; a++) {
            // This assumes vhs prevalence samples target all mature individuals 
            // i.e. maturity = survey vulnerability
            vhsv_pred_numbers(y, 2*a) = disease_samp_vul(a) * N_y_a(y, a) * vhs_immune(y, a);
            vhsv_pred_numbers(y, 2*a+1) = disease_samp_vul(a) * N_y_a(y, a) * vhs_susceptible(y, a);
        }

        vhsv_pred.row(y) = vhsv_pred_numbers.row(y) / vhsv_pred_numbers.row(y).sum();

    }

    dummy_matrix = vhsv_pred;

    // --------------------- FORECAST QUANTITIES ------------------------ //

    vector<Type> N_a_forecast(nage);                   // numbers-at-age forecast
    N_a_forecast.setZero();

    vector<Type> Ntilde_a_forecast(nage);              // mature numbers-at-age forecast
    Ntilde_a_forecast.setZero();

    vector<Type> Ntilde_agecomp_forecast(nage);         // mature numbers agecomp forecast
    Ntilde_agecomp_forecast.setZero();

    vector<Type> waa_forecast(nage);                    // weight-at-age forecast
    waa_forecast.setZero();
    
    vector<Type> winter_survival_forecast(nage);        // winter survival for forecast
    winter_survival_forecast.setZero();
    
    vector<Type> mean_disease_cov(n_covs);              // disease covariate for forecast
    mean_disease_cov.setZero();
    
    Type Btilde_forecast = 0.0;                         // mature biomass forecast

    vector<Type> Btilde_a_forecast(nage);               // mature biomass-at-age forecast
    Btilde_a_forecast.setZero();

    vector<Type> Btilde_agecomp_forecast(nage);         // mature biomass agecomp forecast
    Btilde_agecomp_forecast.setZero();

    Type mdm_forecast = 0.0;                            // mile-days milt forecast

    // forecast disease covariate prevalence
    for(int b = 0; b < n_covs; b++){
        for(int y = nyr - disease_cov_average_years; y < nyr; y++) {
            mean_disease_cov(b) += disease_covs_calc(y, b) / disease_cov_average_years;
        }    
    }
    
    // set up winter survival forecast vector
    for (int a = 0; a < nage; a++) {

        if(a < nage - 1) {
            winter_survival_forecast(a) = exp(-0.5*Z_0_8);
        } else if (a == nage - 1) { 
            winter_survival_forecast(a) = exp(-0.5*Z_9);            
        }

        for (int b = 0; b < n_covs; b++) {
            winter_survival_forecast(a) *= exp(-beta_mortality(b)*mean_disease_cov(b)*mort_age_impact(a, b));
        }
        
    }

    // forecast recruitment
    // take a mean (in log space) of age-2 fish for recruitment forecast
    // enables accounting for age-2 fish harvested in food/bait fishery
    Type mean_log_rec = 0.0;
    for(int y = nyr - recruitment_average_years; y < nyr; y++) {
        mean_log_rec += log(N_y_a(y-1, 2)) / recruitment_average_years;
        // mean_log_rec += log(N_y_a(y, 3)) / recruitment_average_years;
    }

    // project naa matrix forward 1 year
    for (int a = 0; a < nage; a++) {
        
        if (a < 3) {
            
            N_a_forecast(a) = 0;   
            
        } else if(a == 3) {
            
            N_a_forecast(a) = exp(mean_log_rec);
            N_a_forecast(a) *= summer_survival(nyr-1, a-1);
            N_a_forecast(a) -= foodbait_catch(nyr-1, a-1);
            N_a_forecast(a) *= winter_survival_forecast(a-1);
            
        } else if (a >= 4) {

            N_a_forecast(a) = N_y_a(nyr-1, a-1) - spring_removals(nyr-1, a-1);
            N_a_forecast(a) *= summer_survival(nyr-1, a-1);
            N_a_forecast(a) -= foodbait_catch(nyr-1, a-1);
            N_a_forecast(a) *= winter_survival_forecast(a-1);       

        } 
        
        if (a == nage-1) {
            // plus group
            Type N_a_forecast_plus_group = 0.0;
            N_a_forecast_plus_group += N_y_a(nyr-1, a) - spring_removals(nyr-1, a);
            N_a_forecast_plus_group *= summer_survival(nyr-1, a);
            N_a_forecast_plus_group -= foodbait_catch(nyr-1, a);
            N_a_forecast_plus_group *= winter_survival_forecast(a);
            N_a_forecast(a) += N_a_forecast_plus_group;
        }
        
    }
    
    // forecast waa
    for (int a = 0; a < nage; a++) {
        waa_forecast(a) = waa.col(a).tail(waa_average_years).sum() / waa_average_years;
    }

    // pre-fishery spawning biomass forecast
    for (int a = 0; a < nage; a++) {
        Ntilde_a_forecast(a) = maturity(a)*N_a_forecast(a);
        Btilde_a_forecast(a) = Ntilde_a_forecast(a)*waa_forecast(a);
        Btilde_forecast += Btilde_a_forecast(a);
    }

    // agecomp forecasts
    Ntilde_agecomp_forecast = Ntilde_a_forecast / Ntilde_a_forecast.sum();
    Btilde_agecomp_forecast = Btilde_a_forecast / Btilde_a_forecast.sum();

    // forecast proportion of population that are female
    Type perc_female_forecast = 0.0;
    for(int y = nyr - perc_female_forecast_years; y < nyr; y++) {
        perc_female_forecast += perc_female(y) / perc_female_forecast_years;
    }

    // milt forecast
    mdm_forecast = Btilde_forecast - expected_spring_harvest; 
    mdm_forecast *= (1 - perc_female_forecast);
    mdm_forecast /= exp(logmdm_c);

    // ------------------------- LIKELIHOODS ---------------------------- //
    
    // ---- L1: purse-seine age-composition ---- //
    
    // Calculate likelihoods from 1980 to 1984
    // increments plus group from age-5 in 1980 to age-9 in 1984 in a similar 
    // fashion to the numbers-at-age calculation 
    
    Type L1 = 0.0;
    vector<Type> L1_years(nyr);
    L1_years.setZero();
    
    for (int y = 0; y < index_1984; y++) {

        if (seine_ess(y) == -9) continue;               // skip years with no fishery
        int plus_group_index = 5 + y;                   // index of plus group 
        Type seine_plus_group = seine_age_comp.row(y).tail(nage-plus_group_index).sum();

        for (int a = 3; a < plus_group_index; a++) {
            // protects observed age classes with 0% from blowing up likelihood
            if (seine_age_comp(y, a) <= 0) continue;
            if (seine_age_comp_est(y, a) == 0) continue;
            L1_years(y) += seine_age_comp(y, a)*log(seine_age_comp_est(y, a)/seine_age_comp(y, a));
        }

        L1_years(y) += seine_plus_group*log(seine_age_comp_est(y,plus_group_index)/seine_plus_group);
        L1 -= seine_ess(y)*L1_years(y);

    } 
    
    // Calculate likelihoods from 1984 to nyr
    
    for (int y = index_1984; y < nyr; y++) {

        if (seine_ess(y) == -9) continue;              // skip years with no fishery

        for (int a = 3; a < nage; a++) {
            if (seine_age_comp(y, a) <= 0) continue;
            if (seine_age_comp_est(y, a) == 0) continue;
            L1_years(y) += seine_age_comp(y, a)*log(seine_age_comp_est(y, a)/seine_age_comp(y, a));
        }

        L1 -= seine_ess(y)*L1_years(y);

    } 
    

    // ---- L2: spawn survey age-composition ---- //
    
    // Calculate likelihoods from 1980 to 1984
    // increments plus group from age-5 in 1980 to age-9 in 1984 in a similar 
    // fashion to the numbers-at-age calculation 
    
    Type L2 = 0.0;
    vector<Type> L2_years(nyr);
    L2_years.setZero();
    
    for (int y = 0; y < index_1984; y++) {

        if (spawn_ess(y) == -9) continue;               // skip years with no survey
        int plus_group_index = 5 + y;                   // index of plus group 
        Type spawn_plus_group = spawn_age_comp.row(y).tail(nage-plus_group_index).sum();

        for (int a = 3; a < plus_group_index; a++) {
            // protects observed age classes with 0% from blowing up likelihood
            if (spawn_age_comp(y, a) <= 0) continue;
            if (spawn_age_comp_est(y, a) == 0) continue;
            L2_years(y) += spawn_age_comp(y, a)*log(spawn_age_comp_est(y, a)/spawn_age_comp(y, a));
        }

        L2_years(y) += spawn_plus_group*log(spawn_age_comp_est(y,plus_group_index)/spawn_plus_group);
        L2 -= spawn_ess(y)*L2_years(y);

    } 
    
    // Calculate likelihoods from 1984 to nyr
    
    for (int y = index_1984; y < nyr; y++) {

        if (spawn_ess(y) == -9) continue;               // skip years with no survey

        for (int a = 3; a < nage; a++) {
            if (spawn_age_comp(y, a) <= 0) continue;
            if (spawn_age_comp_est(y, a) == 0) continue;
            L2_years(y) += spawn_age_comp(y, a)*log(spawn_age_comp_est(y, a)/spawn_age_comp(y, a));
        }

        L2 -= spawn_ess(y)*L2_years(y);

    } 

    // ---- L3: Number of eggs deposited ---- //

    Type L3 = 0.0;
    vector<Type> L3_total_var(nyr);
    L3_total_var.setZero();


    // number of years in egg index
    int n_egg = 0;

    for (int y = 0; y < nyr; y++) {
        if (egg(y) == -9) continue;                     // skip years with no egg index
        L3_total_var(y) = pow(egg_se(y), 2) + pow(egg_add, 2);
        L3 += log(pow(L3_total_var(y), .5));
        L3 += (pow(log(Ehat_y(y)/egg(y)), 2) / (2*L3_total_var(y)));
        n_egg += 1; 
    }

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
    L5_total_var.setZero();

    // number of years in PWSSC hydroacoustic survey
    int n_pwssc_hydro = 0;

    for (int y = 0; y < nyr; y++) {
        if (pwssc_hydro_se(y) == -9) continue;          // skip years when no PWSSC hydro survey
        L5_total_var(y) = pow(pwssc_hydro_se(y), 2) + pow(pwssc_hydro_add_var, 2);
        L5 += log(pow(L5_total_var(y), .5));
        L5 += (pow(log(Hhat_pwssc_y(y)/pwssc_hydro(y)), 2) / (2*L5_total_var(y)));
        n_pwssc_hydro += 1;
    }

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
    L7_var.setZero();

    for (int y = 0; y < nyr; y++) {
        L7_var(y) = Jhat_y(y) + (pow(Jhat_y(y), 2)/juvenile_overdispersion);
        if (juvenile_survey(y) == -9) continue;
        L7 -= dnbinom2(Type(juvenile_survey(y)), Jhat_y(y), L7_var(y), true);
    }

    // ---- L8: Seroprevalence survey ---- //

    Type L8 = 0.0;

    for (int y = vhsv_obs_start; y < nyr; y++) {
        for (int a = 0; a < nage; a++) {

            if (vhsv_obs(y, 2*a) <= 0) continue;
            if (vhsv_pred(y, 2*a) == 0) continue;
            
            Type x_bin = vhsv_obs(y, 2*a) * vhsv_ess(y);
            Type n_bin = (vhsv_obs(y, 2*a) + vhsv_obs(y, 2*a+1)) * vhsv_ess(y);
            Type p_bin = vhsv_pred(y, 2*a) / (vhsv_pred(y, 2*a) + vhsv_pred(y, 2*a+1));

            // dummy_matrix(y, a) = dbinom(x_bin, n_bin, p_bin, true);

            L8 -= dbinom(x_bin, n_bin, p_bin, true);

        }
    }

    // matrix<Type> dummy_matrix(nyr, 2*nage);               // dummy matrix
    // dummy_matrix.setZero();
    // dummy_matrix = vhsv_pred;

    // ---- build objective function ---- //

    Type negLogLik = 0.0;        
    Type priors = 0.0;
    
    // add likelihood components
    negLogLik += L1;    // seine age comp 
    negLogLik += L2;    // spawn age comp 
    negLogLik += L3;    // egg dep 
    negLogLik += L4;    // ADFG hydro
    negLogLik += L5;    // PWSSC hydro
    negLogLik += L6;    // milt index
    negLogLik += L7;    // juvenile schools
    negLogLik += L8;    // seroprevalence
    
    // add penalties
    negLogLik += 1000*naa_pen;                  // penalizes negative numbers-at-age
    negLogLik += 1000*ntilde_pen;               // penalizes negetive post fishery naa
    negLogLik += 1000*winter_surv_pen;          // penalizes >1 winter survival
    negLogLik += 1000*summer_surv_pen;          // penalizes >1 summer survival

    // add priors
    priors -= dbeta(mat_age3, Type(9.0), Type(11.0), true);
    priors -= dbeta(mat_age4, Type(18.0), Type(2.0), true);
    priors -= dnorm(milt_add_var, Type(0.33), Type(0.10), true);
    priors -= dnorm(adfg_hydro_add_var, Type(0.30), Type(0.08), true);
    priors -= dnorm(pwssc_hydro_add_var, Type(0.32), Type(0.08), true);
    for (int y = 0; y < nyr-vhsv_est_start; y++) {
        priors -= dbeta(vhs_inf_prob(y), Type(1.5), Type(5.0), true);
        priors -= dbeta(vhs_rec_prob(y), Type(4.0), Type(6.0), true);
    }

    negLogLik += priors;

    // ------------------------------------------------------------------ //
    //                         Report Section                             //
    // ------------------------------------------------------------------ //

    // ---- survey quantities ---- //

    REPORT(seine_age_comp_est);
    REPORT(spawn_age_comp_est);
    REPORT(Ehat_y);
    REPORT(Hhat_adfg_y);
    REPORT(Hhat_pwssc_y);
    REPORT(That_y);
    REPORT(Jhat_y);
    REPORT(vhsv_pred);
    
    // ---- population dynamics ---- //
    
    REPORT(age_0);
    REPORT(N_y_a);
    REPORT(maturity);
    REPORT(Ntilde_y_a);
    REPORT(B_y);
    REPORT(Btilde_y);
    REPORT(Btilde_post_y);
    REPORT(spawn_removals);
    REPORT(summer_survival);
    REPORT(winter_survival);
    REPORT(incidence_sp);
    REPORT(fatalities_sp);
    REPORT(seroprev_sp);
    
    // ---- objective function quantities ---- //
    
    // likelihood components
    REPORT(L1);
    REPORT(L2);
    REPORT(L3);
    REPORT(L3_total_var);
    REPORT(L4);
    REPORT(L5);
    REPORT(L5_total_var);
    REPORT(L6);
    REPORT(L7);
    REPORT(L7_var);
    REPORT(L8);
    REPORT(negLogLik);
    
    // penalties
    REPORT(penCount);
    REPORT(naa_pen);
    REPORT(ntilde_pen);
    REPORT(winter_surv_pen);
    REPORT(summer_surv_pen);
    
    // priors
    REPORT(priors);
    
    // ---- forecast quantities ---- //
    REPORT(winter_survival_forecast);
    REPORT(mean_log_rec);
    REPORT(mean_disease_cov);
    REPORT(N_a_forecast);
    REPORT(Ntilde_a_forecast);
    REPORT(Ntilde_agecomp_forecast);
    REPORT(waa_forecast);
    REPORT(Btilde_forecast); 
    REPORT(Btilde_a_forecast); 
    REPORT(Btilde_agecomp_forecast); 
    REPORT(mdm_forecast);    
    
    // ---- dummy objects for debugging ---- //
    REPORT(dummy);
    REPORT(dummy_vector);
    REPORT(dummy_matrix);

    return negLogLik;
  
}
