function [VarianceDecomposition,ShockNames] =...
    unconditional_variance_decomposition(VariableName,oo_,M_,options_)

% INPUTS
%   VariableName          Name of the variable whose variance decomposition
%                         will be saved
%
%
% OUTPUTS 
%    VarianceDecomposition
%    Shock Names

VarianceDecompositionFull = oo_.gamma_y{options_.ar+2};
RowVariableName           = loc(char(M_.endo_names),VariableName);
VarianceDecomposition     = VarianceDecompositionFull(RowVariableName,:);
ShockNames_aux            = char(M_.exo_names); 
ShockNames                = ShockNames_aux(:,5:end);

end