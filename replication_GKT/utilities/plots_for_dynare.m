function plots_for_dynare(model_names_orig,Impulse,Response,n_col,n_row,marks,horizon,name_fig)

% This function creates a graph that combines the impulse responses of
% different models and/or different shocks for teh same model, generates by
% Dynare. If the responses come from different models, it is assumed that
% all the models have the same names for variables and shocks in Dynare.

latex      = 1;
model_names= model_names_orig(:,1)';
v_do_cusum = cell2mat(Response(:,3)');
v_adj      = [cell2mat(Response(:,4)') 1]*100;
v_do_cusum = [v_do_cusum 	0];
v_div_ss   = [zeros(1,11) 	0];

legend_names            = model_names_orig(:,2)';
shocks                  = Impulse(:,1)';
shocks_latex            = Impulse(:,2)';
v_sel_orig              = Response(:,1)';
v_name                  = Response(:,2)';

for i=1:length(model_names)
    name_mat=[model_names{i} '.mat'];
    load([model_names{i} '.mat'])
    eval(['data.M_.' model_names{i} '=M_;'])
    eval(['data.oo_.' model_names{i} '=oo_;'])
end

for j=1:size(shocks,2)
    u_sel=shocks(j);
    if size(model_names,2)>1
        for i=2:size(model_names,2)
            u_sel=[u_sel shocks(j)];
        end
    end
    v_sel  = [v_sel_orig shocks(j)];
    u_name = shocks_latex(j);
    v_name = [v_name shocks_latex(j)];
    
    tt=horizon;
    ir_all=nan(tt,length(v_sel));
    
    figure('Name',char(u_sel(1)),'units','normalized','outerposition',[0 0 1 1])
    set(gcf,'WindowStyle','docked');
    
    for jj=1:length(v_sel)
        subplot(n_row,n_col,jj);
        va = v_sel{jj};
        for ii=1:length(model_names)
            eval(['M_=data.M_.' model_names{ii} ';'])
            eval(['oo_=data.oo_.' model_names{ii} ';'])
            sh = u_sel{ii};
            if isempty(loc(char(M_.endo_names),va))==0
                if v_div_ss(jj)
                    div_ss=abs(oo_.steady_state(loc(char(M_.endo_names),va)));
                else
                    div_ss=1;
                end
                if isfield(oo_.irfs,[va '_eps_' sh])
                    if v_do_cusum(jj)==0
                        ir = getfield(oo_.irfs,[va '_eps_' sh])...
                            *v_adj(jj)/div_ss;
                    elseif v_do_cusum(jj)==1
                        ir = cumsum(getfield(oo_.irfs,[va '_eps_' sh]))...
                            *v_adj(jj)/div_ss;
                    end
                else
                    ir = zeros(1,horizon);
                end
                ir_all(:,ii)=ir(1:tt);
            else
                ir_all(:,ii)=0;
            end
        end
        for ii=1:length(model_names)
            hold on;
            plot(1:tt,ir_all(:,ii),[marks{ii}],'LineWidth',2);
        end
        grid on;
        if latex==1
            va_name=v_name{jj};
            sh_name=u_name{1};
        else
            va_name=va;
            sh_name=sh;
        end
        title(['$' sh_name ' \Rightarrow ' va_name '$'],'interpreter','latex','FontSize',12);
        set(gca,'FontSize',8,'FontName','Times')
        xlim([1 tt]);
        hold off;
    end
    
    h=legend(legend_names);
    set(h,'Position',[0.4 0 0.23 0.03],...
        'Orientation','Horizontal',...
        'Fontsize',8,...
        'Interpreter','Latex');
    legend('boxoff')
    saveplot([name_fig shocks{j}],'graphs of the paper')
end