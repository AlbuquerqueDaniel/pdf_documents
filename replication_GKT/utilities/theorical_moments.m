function [StdVariable,CoefOfAutoPIB,CoefOfAuto1,CoefOfAuto2,CoefOfAuto3,CoefOfAuto4,CoefOfAuto5]...
    = theorical_moments(VariableName,oo_,M_)

%
% INPUTS
%   Variable Name              
%   Order of autocorrelation    
%   oo_
%   M_
%
% OUTPUTS 
%    Std Variable
%    Coef.auto(PIB,X) 
%    Coef.auto1
%    Coef.auto2
%    Coef.auto3
%    Coef.auto4
%    Coef.auto5
%

PosVariable  = loc(char(M_.endo_names),VariableName);
PosPIB       = loc(char(M_.endo_names),'gam_Y_obs');

% Autocorrelation
for i=1:5
    CoefOfAuto  = oo_.autocorr{i}(PosVariable,PosVariable);
    eval(['CoefOfAuto' num2str(i) '= CoefOfAuto;'])
end

%Standar Deviation
StdVariable = sqrt(oo_.var(PosVariable,PosVariable));
StdPIB      = sqrt(oo_.var(PosPIB,PosPIB));

%Correlation between PIB and VariableName
CoefOfAutoPIB = oo_.var(PosPIB,PosVariable) / (StdPIB*StdVariable);
end