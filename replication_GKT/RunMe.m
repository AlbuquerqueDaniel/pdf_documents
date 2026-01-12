% This file replicates the results from the paper "Search Frictions and the 
% Business Cycle in a Small Open Economy DSGE Model", by Juan Guerra-Salas, 
% Markus Kirchner and Rodrigo Tranamil.

tic

clear
close all
clc

addpath([pwd '\utilities'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Making graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('---------------------------------------------------------\n')
fprintf('Stochastic solution and simulation for IRFs...\n')
fprintf('---------------------------------------------------------\n')
pause(3)

%%% Simulations (all estimated params at posterior mean, except shocks)
% Search model
dynare Search_irf console nostrict nolog
rmdir Search_irf s
%rmdir +Search_irf s
% Standard model
dynare Standard_irf console nostrict nolog
rmdir Standard_irf s
%rmdir +Standard_irf s
% Standard model with smaller capital friction
dynare Standard_GAMA_irf console nostrict nolog
rmdir Standard_GAMA_irf s
%rmdir +Standard_GAMA_irf s

clc

%%% Figure 1 and 2
fprintf('---------------------------------------------------------\n')
fprintf('Producing figure 1 and 2...\n')
fprintf('---------------------------------------------------------\n')

%               Name                  Latex Name
model_names = {'Search_irf_results',  'Search'; 
               'Standard_irf_results','Standard'};

%               Name     Latex Name
Impulse     = {'Rstar',  'R^{\ast}';
               'eR',     'e^R'};       
        
%               Name       Latex Name   Cusum    Factor
Response    = {'gam_YNM'   ,'Y^{NM}'    ,1       ,1;
               'gam_C'     ,'C'         ,1       ,1;
               'gam_I'     ,'I'         ,1       ,1;
               'rer'       ,'RER'       ,0       ,1;               
               'h'         ,'h'         ,0       ,1;
               'rho'       ,'\rho'      ,0       ,1; 
               'pi'        ,'4\pi'      ,0       ,4;
               'gam_W'     ,'4\pi^{W}'  ,0       ,4;               
               'u'         ,'u'         ,0       ,1;
               'v'         ,'v'         ,0       ,1;
               'thetaL'    ,'\theta'    ,0       ,1}; 
n_col    = 4;
n_row    = 3;
marks    = {'b','-.r'};
horizon  = 20;
name_fig = 'irf_';
plots_for_dynare(model_names,Impulse,Response,n_col,n_row,marks,...
    horizon,name_fig)
pause(3)
close all
clc

%%% Figure 3
fprintf('---------------------------------------------------------\n')
fprintf('Producing figure 3...\n')
fprintf('---------------------------------------------------------\n')

%               Name                        Latex Name
model_names = {'Standard_irf_results',      'Standard'; 
               'Standard_GAMA_irf_results', 'Standard Smaller $\gamma$'};

%               Name     Latex Name
Impulse     = {'Rstar',  'R^{\ast}';
               'eR',     'e^R'};       
        
%               Name       Latex Name   Cusum    Factor
Response    = {'gam_YNM'   ,'Y^{NM}'    ,1       ,1;
               'gam_C'     ,'C'         ,1       ,1;
               'gam_I'     ,'I'         ,1       ,1;
               'h'         ,'h'         ,0       ,1;
               'pi'        ,'4\pi'      ,0       ,4};  
n_col    = 2;
n_row    = 3;
marks    = {'-.r','--g'};
horizon  = 20;
name_fig = 'irf_GAMA_';
plots_for_dynare(model_names,Impulse,Response,n_col,n_row,marks,...
    horizon,name_fig)
pause(3)

close all
clc

%%% Figure 4 and 5
fprintf('---------------------------------------------------------\n')
fprintf('Producing figure 4 and 5...\n')
fprintf('---------------------------------------------------------\n')

model_names      = {'Standard','Search','BVAR1','BVAR2','BVAR3'};
vars_names       = {'gam_YNM_obs','pi_obs','R_obs','rer_obs','hn_obs','gam_W_obs'};
vars_names_latex = {'\Delta \log Y^{NM}','\pi','R','rer','h\times n','\Delta \log W'};
recursive_forecast_for_plot(model_names,vars_names)
plot_recursive_forecasts(model_names,vars_names,vars_names_latex)
pause(3)

close all
clc

fprintf('---------------------------------------------------------\n')
fprintf('Figures done.\n')
fprintf('---------------------------------------------------------\n')
pause(3)

clc

fprintf('---------------------------------------------------------\n')
fprintf('Making tables now...\n')
fprintf('---------------------------------------------------------\n')
pause(3)

clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Making tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation of each model (all estimated params at posterios mean)
% Search model
dynare Search console nostrict nolog
Search.oo_ = oo_;
Search.M_ = M_;
Search.options_ = options_;
rmdir Search s
%rmdir +Search s
% Standard model
dynare Standard console nostrict nolog
Standard.oo_ = oo_;
Standard.M_ = M_;
Standard.options_ = options_;
rmdir Standard s
%rmdir +Standard s

clc

fprintf('---------------------------------------------------------\n')
fprintf('Table 3: Second Moments.\n')
fprintf('---------------------------------------------------------\n')

Variables  = {'gam_YNM_obs','gam_C_obs','gam_I_obs','pi_obs','R_obs',...
    'rer_obs','xi_obs','gam_W_obs','hn_obs','u_obs'};
VariableNames  = {'Non-mining GDP','Consumption','Investment','Inflation',...
    'MPR','Real exch. rate','EMBIG Chile','Real wage','Total hours worked',...
    'Unemployment rate'};

% Table 3 Search model
Std = [];
AC1 = [];

for i=1:length(Variables)
    [StdVariable,~,CoefOfAuto1,~,~,~,~] =...
        theorical_moments(Variables{i},Search.oo_,Search.M_);
    Std =  [Std StdVariable];
    AC1 =  [AC1 CoefOfAuto1];
end

title   = 'Second moments Search model:';
headers = char('','Standard Deviation','AC Order 1');
labels  = char(VariableNames);
values = [Std',AC1'];
dyntable2(title,headers,labels,values,0,5,2)

% Table 3 Standard model
Std = [];
AC1 = [];

for i=1:length(Variables)
    [StdVariable,~,CoefOfAuto1,~,~,~,~] =...
        theorical_moments(Variables{i},Standard.oo_,Standard.M_);
    Std =  [Std StdVariable];
    AC1 =  [AC1 CoefOfAuto1];
end

title   = 'Second moments Standard model:';
headers = char('','Standard Deviation','AC Order 1');
labels  = char(VariableNames);
values = [Std',AC1'];
dyntable2(title,headers,labels,values,0,5,2)

fprintf(' \n')
fprintf('---------------------------------------------------------\n')
fprintf('Table 5: Variance Decomposition.\n')
fprintf('---------------------------------------------------------\n')

Variables_FEVD  = {'yNM','c','i','pi','R',...
    'rer','xi_obs','w','hn_obs','u_obs'};

% Table 5 Search model
VarianceDecomposition = [];
for i=1:length(Variables_FEVD)
    [VarianceDecomp,~] =...
        unconditional_variance_decomposition(Variables_FEVD{i},Search.oo_,Search.M_,Search.options_);
    VarianceDecomposition = [VarianceDecomposition; VarianceDecomp];
end

Tech    = VarianceDecomposition(:,3)+VarianceDecomposition(:,5)+VarianceDecomposition(:,6);
Pref    = VarianceDecomposition(:,1)+VarianceDecomposition(:,2);
MPRate  = VarianceDecomposition(:,15);
GovCons = VarianceDecomposition(:,8);
Foreign = VarianceDecomposition(:,9)+VarianceDecomposition(:,10)+VarianceDecomposition(:,13)+VarianceDecomposition(:,14);
Risk    = VarianceDecomposition(:,11)+VarianceDecomposition(:,12);
CoProd  = VarianceDecomposition(:,7);
ExoSep  = VarianceDecomposition(:,4);
VarDecom = [Tech Pref MPRate GovCons Risk Foreign CoProd ExoSep];

title   = 'Variance Decomposition (in percentage) Search Model';
headers = char('','Tech.','Pref.','MP Rate','Gov. Cons.','Risk','Foreign','Co. Prod.','Exo. Sep.');
labels  = char(VariableNames);
values = round(100*VarDecom);
dyntable2(title,headers,labels,values,0,0,0)

% Table 5 Standard model
VarianceDecomposition = [];
for i=1:length(Variables_FEVD)
    [VarianceDecomp,ShockNames] =...
        unconditional_variance_decomposition(Variables_FEVD{i},Standard.oo_,Standard.M_,Standard.options_);
    VarianceDecomposition = [VarianceDecomposition; VarianceDecomp];
end

Tech    = VarianceDecomposition(:,3)+VarianceDecomposition(:,4)+VarianceDecomposition(:,5);
Pref    = VarianceDecomposition(:,1)+VarianceDecomposition(:,2);
MPRate  = VarianceDecomposition(:,8);
GovCons = VarianceDecomposition(:,14);
Foreign = VarianceDecomposition(:,10)+VarianceDecomposition(:,11)+VarianceDecomposition(:,12)+VarianceDecomposition(:,13);
Risk    = VarianceDecomposition(:,6)+VarianceDecomposition(:,7);
CoProd  = VarianceDecomposition(:,9);
VarDecom = [Tech Pref MPRate GovCons Risk Foreign CoProd];

title   = 'Variance Decomposition (in percentage) Standard Model';
headers = char('','Tech.','Pref.','MP Rate','Gov. Cons.','Risk','Foreign','Co. Prod.');
labels  = char(VariableNames);
values = round(100*VarDecom);
dyntable2(title,headers,labels,values,0,0,0)

clearvars
disp(' ')
toc;