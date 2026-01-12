clear;
clc;
close all;

DIR='figs1-2';

%                mat file                       legend name
model_names = {	'Search_results',               'Search'; %Paper and appendix
                'Standar_results',              'Standard'; %Paper and appendix
		        %'GKT2019_results_base_stochsim_corr_priormeanshks_slackness',    'Search'; %Paper and appendix
                %'CalvoW_stochsimul_postmean_shockspriormean_results',             'Standard'; %Paper and appendix
		        %'CalvoW_stochsim_postmean_gammasearch_shockspriormean_results',   'Standard low GAMA'; %Paper and appendix
                %'GKT2019_results_lowendo_stochsim_priormeansshks_slackness',     'Search Lower Endo. Sep.'; %Appendix
              };

%                Name        Latex name        shock scale
Impulse     = {	%'varrho'     ,'\varrho'          , 1;
                %'kappa'      ,'\kappa'           , 1;
                %'varomega'   ,'\omega'           , 1;
                %'g'          ,'g'                , 1;
                %'z'          ,'z'                , 1;
                %'a'          ,'a'                , 1;
                %'zetao'      ,'\zeta^o'          , 1;
                %'zetau'      ,'\zeta^u'          , 1;              
                'eR'         ,'e^R'              , 1;
                %'yCo'        ,'y^{Co}'           , 1;
                'Rstar'      ,'R^{\ast}'         , 1;
                %'pistar'     ,'\pi^{\ast}'       , 1;
                %'pCostar'    ,'p^{{Co}^\ast}'    , 1;
                %'ystar'      ,'y^{\ast}'         , 1;
                %n'm'          ,'m'                , 1;
                %'rhox'       ,'\rho^x'           , 1;
                %'rho'        ,'\rho'             , 1;
              };       
        
%                Name         Latex Name   Cusum    Factor    Div_ss
Response    = {'gam_YNM'      ,'Y^{NM}'    ,1       ,1        ,0;
               'gam_C'        ,'C'         ,1       ,1        ,0;
               'gam_I'        ,'I'         ,1       ,1        ,0;
               'rer'          ,'RER'       ,0       ,1        ,0;
               
               'h'            ,'h'         ,0       ,1        ,0;
               'rho'          ,'\rho'      ,0       ,1        ,0; 
               'pi'           ,'4\pi'      ,0       ,4        ,0;
               'gam_W'        ,'4\pi^{W}'  ,0       ,4        ,0;
                 
               'u'            ,'u'         ,0       ,1        ,0;
               'v'            ,'v'         ,0       ,1        ,0;
               'thetaL'        ,'\theta'   ,0       ,1        ,0;         
              };
          
plot_options.n_col      = 4;
plot_options.n_row      = 3;

plot_options.marks      = {'b','-.r','--g','-.k'};
%plot_options.marks      = {'-.r','--g'};

plot_options.model_names= model_names(:,1)';
legend_names            = model_names(:,2)';

shocks                  = Impulse(:,1)';
shocks_latex            = Impulse(:,2)';
shocks_scale            = Impulse(:,3)';

v_sel                   = Response(:,1)';
v_name                  = Response(:,2)';
plot_options.v_do_cusum = cell2mat(Response(:,3)');
plot_options.v_adj      = cell2mat(Response(:,4)');
plot_options.v_div_ss   = cell2mat(Response(:,5)');


plot_options.v_do_cusum = [plot_options.v_do_cusum 	0];  	% shock doesn't sum over
plot_options.v_div_ss   = [plot_options.v_div_ss 	0];    	% shock doesn't divided by its ss
plot_options.v_adj_orig = [plot_options.v_adj  		1]*100;	% adjustment for shock

plot_options.horizon    = 20;
plot_options.grid       = 1;
plot_options.latex      = 1;

mkdir(DIR)

for ii=1:length(plot_options.model_names)
    name_mat=[plot_options.model_names{ii} '.mat'];
    load([plot_options.model_names{ii} '.mat'])
    eval(['data.M_.' plot_options.model_names{ii} '=M_;'])
    eval(['data.oo_.' plot_options.model_names{ii} '=oo_;'])
end

for j=1:size(shocks,2)
    plot_options.u_sel=shocks(j);
    if size(plot_options.model_names,2)>1
        for i=2:size(plot_options.model_names,2)
            plot_options.u_sel=[plot_options.u_sel shocks(j)];
        end
    end
    plot_options.v_sel  = [v_sel shocks(j)];
    plot_options.u_name = shocks_latex(j);
    plot_options.v_name = [ v_name shocks_latex(j)];
    plot_options.v_adj  = shocks_scale{j} * plot_options.v_adj_orig;
    
    plots_for_dynare(plot_options,data)
    
    h=legend(legend_names);
    set(h,'Position',[0.4 0 0.23 0.03],...
        'Orientation','Horizontal',...
        'Fontsize',8,...
        'Interpreter','Latex');
    legend('boxoff')

    fig_ratio=16/9;
    trim_hor=0.6;
    trim_top=0.05;
    saveplot(['irf_' shocks{j}],fig_ratio,trim_hor,trim_top,DIR)
end