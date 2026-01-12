clear;
clc;

% Define model name(s) to work  (use dir *.mat)
ModelNames  = {'Search','Standar'};

% Define variable(s) to work
VariableNames  = {'gam_YNM_obs','gam_C_obs','gam_I_obs','pi_obs','R_obs',...
    'rer_obs','xi_obs','gam_W_obs','hn_obs','u_obs'};
UseAllVarobs    = 0; % 1 = all varobs instead of VariableNames


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go through models
for j=1:length(ModelNames)
    disp('---------------------------')
    disp(['Summary for ' ModelNames{j} ' Model'])
    disp('---------------------------')
    
    %%%%%%% Load results of model j %%%%%%%
    load([ModelNames{j} '_results.mat'],'oo_','M_','options_')
    
    %%%%%%% Find varobs of the model j %%%%%%%
    if UseAllVarobs == 1
        VariableNames = cellstr(options_.varobs);
    end
    
    %%%%%%% Table 2: Marginal Data Density %%%%%%%
    
    %%%%%%% Table 3: Second Moments %%%%%%%
    Std = [];
    CoefOfAutoPIB = [];
    Autocorrelation1 = [];
    Autocorrelation2 = [];
    Autocorrelation3 = [];
    Autocorrelation4 = [];
    Autocorrelation5 = [];
    
    for i=1:length(VariableNames)
        [StdVariable,CoefAutoPIB,CoefOfAuto1,CoefOfAuto2,CoefOfAuto3,CoefOfAuto4,CoefOfAuto5] =...
            theorical_moments(VariableNames{i},oo_,M_);
        Std              =  [Std StdVariable];
        CoefOfAutoPIB    =  [CoefOfAutoPIB CoefAutoPIB];
        Autocorrelation1 =  [Autocorrelation1 CoefOfAuto1];
        Autocorrelation2 =  [Autocorrelation2 CoefOfAuto2];
        Autocorrelation3 =  [Autocorrelation3 CoefOfAuto3];
        Autocorrelation4 =  [Autocorrelation4 CoefOfAuto4];
        Autocorrelation5 =  [Autocorrelation5 CoefOfAuto5];
    end
    
    title   = ['Table 3: Second Moments - ' ModelNames{j} ' Model'];
    headers = char('','Standard Deviation','Correlation with GDP','Autocorrel(1)',...
        'Autocorrel(2)','Autocorrel(3)','Autocorrel(4)','Autocorrel(5)');
    labels  = char(VariableNames);
    values = [Std',CoefOfAutoPIB',Autocorrelation1',...
        Autocorrelation2',Autocorrelation3',...
        Autocorrelation4',Autocorrelation5'];
    dyntable(title,headers,labels,values,0,5,2)
    
    %%%%%%% Table 4: Estimated Parameters %%%%%%%
    [ParameterNames,ParameterPriors,ParameterPostMode,ParameterPostSD] =...
        estimated_parameters(oo_);
    
    title   = ['Table 4: Estimated Parameters - ' ModelNames{j} ' Model'];
    headers = char('','Prior Mean','Posterior Mean','Posterior s.d.');
    labels  = char(ParameterNames);
    values  = [ParameterPriors ParameterPostMode ParameterPostSD];
    dyntable(title,headers,labels,values,0,0,4)
    
    %%%%%%% Table 5: Variance Decomposition %%%%%%%
    VarianceDecomposition = [];
    for i=1:length(VariableNames)
        [VarianceDecomp,ShockNames] =...
            unconditional_variance_decomposition(VariableNames{i},oo_,M_,options_);
        VarianceDecomposition = [VarianceDecomposition; VarianceDecomp];
    end
    
    title   = ['Table 5: Variance Decomposition (in percentage) - ' ModelNames{j} ' Model'];
    headers = char('',ShockNames);
    labels  = char(VariableNames);
    values = 100*VarianceDecomposition;
    dyntable(title,headers,labels,values,0,0,1)
    disp(' ')
    disp(' ')
end