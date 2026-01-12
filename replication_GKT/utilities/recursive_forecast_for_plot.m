function recursive_forecast_for_plot(model_names,vars_names)

% Data to be used
%[~, ~, raw] = xlsread([pwd '\estim_data.xlsx'],'estim_data','A2:AH60');
[~, ~, raw] = xlsread('.\estim_data.xlsx','estim_data','A2:AH60');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
r = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw);
raw(r) = {NaN};
data = reshape([raw{:}],size(raw));

actual.dates = data(:,1);
actual.gam_YNM = data(:,2);
actual.pi    = data(:,13);
actual.R     = data(:,15);
actual.rer   = data(:,17);
actual.hn    = data(:,31);
actual.gam_W = data(:,29);

% Setting of forecast
param_fct.k = 22;                           % Size of first recursive sample
param_fct.h = 10;                           % # of forecast periods
param_fct.s = 59;                           % Size of total sample
param_fct.T = param_fct.s+param_fct.h;      % Size of total sample + # of forecast periods

% Matrix to be used
mat_smoothed=NaN(param_fct.s,param_fct.s-param_fct.k+1);
mat_forecast=NaN(param_fct.s,param_fct.s-param_fct.k+1);
mat_total=NaN(param_fct.T,param_fct.s-param_fct.k+1);
mat_error=NaN(param_fct.s,param_fct.s-param_fct.k+1);
mat_error_aux=NaN(param_fct.h,param_fct.s-param_fct.k+1);
mat_RMSE=NaN(param_fct.h,1);

for m=1:length(model_names)
    for j=1:length(vars_names)
        var_name=vars_names{j};
        var_actual=eval(strcat('actual.',var_name(1:end-4)));
        do_loop=1;
        i=1;
        
        while do_loop
            model.num  = num2str(i+param_fct.k-1);
            model.ext  = '.mat';
            model.name = strcat([pwd '\recursive forecasts\estimation_models\' model_names{m}],'_',model.num,model.ext);
            if exist(model.name,'file')==2
                %Smoothed variables (always Search model)
                model_recur=load(strcat([pwd '\recursive forecasts\estimation_models\' 'Search'],'_',model.num,model.ext));
                var_ss=model_recur.oo_.steady_state(loc(char(model_recur.M_.endo_names),var_name));
                var_smoothed_aux=eval(['model_recur.oo_.SmoothedVariables.' var_name]);
                var_smoothed=[var_smoothed_aux; NaN(param_fct.T-length(var_smoothed_aux),1)];
                
                %Forecasted variables
                if strcmp(model_names{m},'Standard') || strcmp(model_names{m},'Search')
                    model_recur=load(model.name);
                    var_forecast=eval(['model_recur.oo_.forecast.Mean.' var_name])-var_ss;
                else
                    model_recur=load(model.name);
                    var_forecast=eval(['model_recur.oo_.bvar.forecast.no_shock.Mean.' var_name])-var_ss;
                end
                var_forecast_aux=[NaN(size(var_smoothed_aux)); var_forecast];
                var_forecast_aux(length(var_smoothed_aux))=var_smoothed(length(var_smoothed_aux)); %Adding starting point
                var_forecast=[var_forecast_aux; NaN(param_fct.s-length(var_forecast_aux),1)];
                
                %Matrix containing smoothed variables for each period
                mat_smoothed(:,i)=var_smoothed(1:param_fct.s);
                %Matrix containing forecasted variables for each period
                mat_forecast(:,i)=var_forecast(1:param_fct.s);
                if i<=size(mat_total,2)
                    %Matrix containing smoothed+forecasted variables for each period
                    mat_total(:,i)=[var_smoothed_aux ; var_forecast(length(var_smoothed_aux)+1:end);NaN(length(mat_total)-length(var_forecast),1)];
                else
                    break
                end
            else
                break
            end
            i=i+1;
        end
        
        % Compute RMSE
        for k=1:i-1
            if k<=param_fct.s-param_fct.k+1
                mat_error(:,k)=var_actual-mat_forecast(:,k);
                mat_error_aux(:,k)=[NaN(param_fct.k-1+k+param_fct.h-param_fct.s,1); mat_error(param_fct.k+k:min(param_fct.k-1+k+param_fct.h,param_fct.s),k)];
            else
                break
            end
        end
        mat_RMSE=sqrt(mean(mat_error_aux.^2,2,'omitnan'));
        std_actual = std(var_actual);

        % Save results in MAT-file
        save(strcat([pwd '\recursive forecasts\mat_for_plot\' model_names{m}],'_',var_name,'.mat'),'mat_smoothed','mat_forecast','mat_total','mat_RMSE','var_actual','std_actual') 
    end   
end
end

