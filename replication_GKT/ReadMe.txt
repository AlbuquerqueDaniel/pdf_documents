Contents of this folder:

1. Main Matlab script RunMe.m.
2. Archive estim_data.xls.
3. Various Dynare model-files (*.mod).
4. Sub-folder "utilities" with Matlab functions.
5. Sub-folder "MH estimation" with various Dynare model-files to perform MH estimation.
6. Sub-folder "recursive forecasts" with various Dynare model-files to perform recursive forecasts.
7. Sub-folder "graphs of the paper" containing the figures of the paper in pdf, eps-format.


************************** Results without estimation ************************************************

Run the file RunMe.m to replicate the results of the paper. Make sure that Matlab and Dynare are installed 
correctly on your computer. For details about the installation process of these software please see 
mathworks.com for Matlab and www.dynare.org for Dynare.

The produced graphs are saved in the folder "graphs of the paper". 
The mapping is as follows:
	irf_Rstar.eps      --> figure 1 (Impulse Responses to a Foreign Interest Rate Shock) 
	irf_eR.eps         --> figure 2 (Impulse Responses to a Domestic Monetary Policy Shock)
	irf_GAMA_eR.eps    --> figure 3 (Impulse Responses of the Standard Model with a Smaller Capital Friction) 
	irf_GAMA_Rstar.eps --> figure 3 (Impulse Responses of the Standard Model with a Smaller Capital Friction) 
	Forecasts.eps      --> figure 4 (Recursive Out-of-Sample Forecasts)
	RMSE.eps           --> figure 5 (Root Mean Squared Errors of Out-of-Sample Forecasts)
 
The inputs to construct the tables of the paper are printed to the screen. 
The mapping is as follows:
	Second moments Search model    --> Column 4 and 7 of table 3
	Second moments Standard model  --> Column 5 and 8 of table 3
	Variance Decomp Search model   --> Column 1-8, 1st panel of table 5
	Variance Decomp Standard model --> Column 1-8, 2nd panel of table 5

The results in the paper were generated with Dynare version 4.6.1 and Matlab version 9.8.0.1396136 (R2020a) Update 3. 
Total runtime (without estimation) on a laptop with an Intel Core i5-1035G1 CPU @ 1.00GHz with 8 GB RAM 
is about 70 seconds.

***************************  Estimation Procedure ****************************************************

The MH estimation of each model can be done by running:
	...\replication\MH estimation\Search.mod
	...\replication\MH estimation\Standard.mod
The estimation procedure uses the posterior mode as the initial point for the MH algorithm. The Dynare 
output of each estimation gives you information about the posterior mean of all estimated parameters and 
the marginal data density of the estimation. Note that estimated parameters could be slightly different 
from the paper because of the random number generation and/or Dynare version you are using. To generate 
the graphs and tables, you need to move each *_results.mat to the main folder and then run "RunMe.m".

The recursive forecasts of each model can be done by running:
	...\replication\recursive forecasts\estimation_models\BVAR1.mod
	...\replication\recursive forecasts\estimation_models\BVAR2.mod
	...\replication\recursive forecasts\estimation_models\BVAR2.mod
	...\replication\recursive forecasts\estimation_models\Search.mod
	...\replication\recursive forecasts\estimation_models\Standard.mod
Each Dynare file generate mat-files that are used by "RunMe.m" to construct the recursive-forecast graphs. 
You don't need to do nothing else, "RunMe.m" knows where these mat-file are.

Total running time is several hours (at least 10 hours).