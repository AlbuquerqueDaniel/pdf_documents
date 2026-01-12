function [ParameterNames,ParameterPriors,ParameterPostMode,ParameterPostSD] =...
    estimated_parameters(oo_)

% INPUTS
%   oo_            
%
% OUTPUTS 
%    Parameter Names
%    Parameter Priors
%    Parameter Post Mode
%    Parameter Post s.d. 

ParameterNames1      = fieldnames(oo_.posterior_mode.parameters);
ParameterValues1     = zeros(size(ParameterNames1,1),1);
ParameterPriorMean1  = oo_.prior.mean(end-size(ParameterNames1,1)+1:end,1);
ParameterPriorSD1    = diag(oo_.posterior.optimization.Variance(end-size(ParameterNames1,1)+1:end,...
                        end-size(ParameterNames1,1)+1:end)).^0.5;
for i=1:size(ParameterNames1,1)
    ParameterValues1(i) = getfield(oo_.posterior_mode.parameters,ParameterNames1{i});
end

ParameterNames2      = fieldnames(oo_.posterior_mode.shocks_std);
ParameterValues2     = zeros(size(ParameterNames2,1),1);
ParameterPriorMean2  = oo_.prior.mean(1:size(ParameterNames2,1),1);
ParameterPriorSD2    = diag(oo_.posterior.optimization.Variance(1:size(ParameterNames2,1),...
                        1:size(ParameterNames2,1))).^0.5;
for i=1:size(ParameterNames2,1)
    ParameterValues2(i) = getfield(oo_.posterior_mode.shocks_std,ParameterNames2{i});
end
 
ParameterNames     = [ParameterNames1;     ParameterNames2];
ParameterPriors    = [ParameterPriorMean1; ParameterPriorMean2];
ParameterPostMode  = [ParameterValues1;    ParameterValues2];
ParameterPostSD    = [ParameterPriorSD1;   ParameterPriorSD2];

end