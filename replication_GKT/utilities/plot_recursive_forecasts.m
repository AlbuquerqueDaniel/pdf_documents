function plot_recursive_forecasts(model_names,vars_names,vars_names_latex)
% Plot_recursive

% Setting of forecast
k = 22;       % Size of first recursive sample
h = 10;       % # of forecast periods
s = 59;       % Size of total sample
T = s+h;      % Size of total sample + # of forecast periods
dates = (2001.50:0.25:2016.00)';  % Dates


for j=1:length(vars_names)
    %Plot figure for Search recursive forecasts
    load(strcat([pwd '\recursive forecasts\mat_for_plot\' 'Search'],'_',vars_names{j}))
    figure(1)
    set(gcf,'WindowStyle','docked')
    subplot(length(vars_names),2,2*j-1)
    box on
    plot(dates,var_actual,'k','LineWidth',1.5);
    hold on
    plot(dates,mat_forecast,'r');
    xlim([2001.5 2016.00])
    ylabel(['$' vars_names_latex{j} '$'],'interpreter','latex','FontName','Times','FontSize',50,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
    if j==1
        title('Search','FontName','Times','FontSize',12);
    end
    set(gca,'FontSize',12)
    set(gcf, 'PaperPosition', [0 0 14 14]);
    set(gcf, 'PaperSize', [14 14]);
    
    %Plot figure for Standard recursive forecast
    load(strcat([pwd '\recursive forecasts\mat_for_plot\' 'Standard'],'_',vars_names{j}))
    subplot(length(vars_names),2,2*j);
    plot(dates,var_actual,'k','LineWidth',1.5);
    hold on
    plot(dates,mat_forecast,'r');
    xlim([2001.5 2016.00])
    if j==1
        title('Standard','FontName','Times','FontSize',12)
    end
    set(gca,'FontSize',12)
    set(gcf, 'PaperPosition', [0 0 14 14]);
    set(gcf, 'PaperSize', [14 14]);
    
    if j==length(vars_names)
        leg=legend('Data','Forecasts');
        set(leg,'Position',[0.329 0.011 0.370 0.064],...
            'Orientation','Horizontal',...
            'FontName','Times',...
            'Fontsize',12);
        saveas(gcf,[pwd '\graphs of the paper\forecasts.eps'],'epsc2');
        saveas(gcf,[pwd '\graphs of the paper\forecasts.pdf'],'pdf');
    end
    
    %Plot RMSE
    figure(2)
    set(gcf,'WindowStyle','docked')
    subplot(ceil(sqrt(length(vars_names))),length(vars_names)/ceil(sqrt(length(vars_names))),j)
    box on
    horizon=1:h;
    load(strcat([pwd '\recursive forecasts\mat_for_plot\' 'BVAR1'],'_',vars_names{j}))
    plot(horizon,mat_RMSE,'-ob');
    hold on
    load(strcat([pwd '\recursive forecasts\mat_for_plot\' 'BVAR2'],'_',vars_names{j}))
    plot(horizon,mat_RMSE,'-*m');
    load(strcat([pwd '\recursive forecasts\mat_for_plot\' 'BVAR3'],'_',vars_names{j}))
    plot(horizon,mat_RMSE,'-xg');
    load(strcat([pwd '\recursive forecasts\mat_for_plot\' 'Search'],'_',vars_names{j}))
    plot(horizon,mat_RMSE,'k','LineWidth',1.5);
    load(strcat([pwd '\recursive forecasts\mat_for_plot\' 'Standard'],'_',vars_names{j}))
    plot(horizon,mat_RMSE,'--r');
    xlim([1 h]);
    title(['$' vars_names_latex{j} '$'],'interpreter','latex','FontName','Times','FontSize',50)
    set(gca,'FontSize',18)
    set(gcf, 'PaperPosition', [0 0 20 20]);
    set(gcf, 'PaperSize', [20 20]);
    if j==length(vars_names)
        leg=legend('BVAR1','BVAR2','BVAR3','Search','Standard'); %,'s.d. Data'
        set(leg,'Position',[0.4 0.01 0.23 0.03],...
            'Orientation','Horizontal',...
            'FontName','Times',...
            'Fontsize',18);
        saveas(gcf,[pwd '\graphs of the paper\RMSE.eps'],'epsc2');
        saveas(gcf,[pwd '\graphs of the paper\RMSE.pdf'],'pdf');
    end
end

end

