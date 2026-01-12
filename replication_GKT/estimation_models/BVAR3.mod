close all;

//BVARs
@#define BVAR1 = 0
@#define BVAR2 = 0
@#define BVAR3 = 1

//Periods
@#define periods = 22:59

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Variables
@#if BVAR1
    var gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs;
@#endif
@#if BVAR2
    var gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs 
        gam_C_obs gam_I_obs g_obs;
@#endif
@#if BVAR3
    var gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs 
        ystar_obs Rstar_obs pCostar_obs pistar_obs xi_obs;
@#endif

//Observed variables
@#if BVAR1
    varobs gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs;
@#endif
@#if BVAR2
    varobs gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs 
        gam_C_obs gam_I_obs g_obs;
@#endif
@#if BVAR3
    varobs gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs 
        ystar_obs Rstar_obs pCostar_obs pistar_obs xi_obs;
@#endif

// Recursive estimation 
@#for n in periods  
    bvar_forecast(datafile=estim_data, xls_sheet=estim_data, first_obs=5,
                  forecast=10, bvar_replic=10000, nobs=@{n-4}, prefilter=0) 4;
    @#if BVAR1 
        save BVAR1_@{n}.mat
    @#endif
    @#if BVAR2 
        save BVAR2_@{n}.mat
    @#endif
    @#if BVAR3 
        save BVAR3_@{n}.mat
    @#endif
    close all
@#endfor