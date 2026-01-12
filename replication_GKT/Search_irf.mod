% NKSOE model with staggered nash wage bargaining and endogenous separations

%--------------------------------------------------------------------------
% 0. Housekeeping (close all graphic windows) and options
%--------------------------------------------------------------------------
close all;

%--------------------------------------------------------------------------
% 1. Preamble
%-------------------------------------------------------------------------
%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%
//Engogenous variables
var lam q c i k theta chitil u m n rho s e thetaL Ebar Estar U Delbar 
Delstar wstar epsi mu Fstar Wstar Wbar w wbar wtil h yC xH xF yH yF pm 
ptilH ptilF fH fF mcH mcF rK v rhon ctil hC imp pi R xi rer piS xHstar pH 
pF DelH DelF dstar tb y yNM pY pYNM;

//Exogenous state variables
var kappa varrho varomega rhox z a yCo g Rstar pistar zetao zetau ystar 
pCostar eR;

//Exogenous innovations
varexo eps_kappa eps_varrho eps_varomega eps_rhox eps_z eps_a eps_yCo eps_g
eps_Rstar eps_pistar eps_zetao eps_zetau eps_ystar eps_pCostar eps_eR eps_aux;

//Definitions (4)
var stb gam_Y gam_C gam_I gam_W pCoyCo gam_YNM;

//Observed variables
var gam_Y_obs gam_YNM_obs gam_C_obs gam_I_obs gam_W_obs R_obs xi_obs pi_obs
rer_obs hn_obs u_obs pCoyCo_obs g_obs Rstar_obs pistar_obs ystar_obs pCostar_obs; 

%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%
//Deep parameters
parameters VARSIGMA GAMA BETTA SIGMA DELTA UPSILON MU ALPHA ETA ETASTAR
ALPHA_W VARTHETA_W THETA_W BBAR PHI OMEGA VARPHI SIGMA_CTIL MU_CTIL MBAR
VARTHETA_H THETA_H EPSILON_H VARTHETA_F THETA_F EPSILON_F OSTAR O CHI PSI 
RHO_R ALPHA_PI ALPHA_Y;

//AR parameters
parameters RHO_kappa RHO_varrho RHO_varomega RHO_rhox RHO_z RHO_a RHO_yCo 
RHO_g RHO_Rstar RHO_pistar RHO_zetao RHO_zetau RHO_ystar RHO_pCostar RHO_eR;

parameters SIG_kappa SIG_varrho SIG_varomega SIG_rhox SIG_z SIG_a SIG_yCo 
SIG_g SIG_Rstar SIG_pistar SIG_zetao SIG_zetau SIG_ystar SIG_pCostar SIG_eR;

//SS parameters
parameters lam_ss q_ss c_ss i_ss k_ss theta_ss chitil_ss u_ss m_ss n_ss 
rho_ss s_ss e_ss thetaL_ss Ebar_ss Estar_ss U_ss Delbar_ss Delstar_ss 
wstar_ss epsi_ss mu_ss Fstar_ss Wstar_ss Wbar_ss w_ss wbar_ss wtil_ss h_ss
yC_ss xH_ss xF_ss yH_ss yF_ss pm_ss ptilH_ss ptilF_ss fH_ss fF_ss mcH_ss 
mcF_ss rK_ss v_ss rhon_ss ctil_ss hC_ss imp_ss pi_ss R_ss xi_ss rer_ss 
piS_ss xHstar_ss pH_ss pF_ss DelH_ss DelF_ss dstar_ss tb_ss y_ss yNM_ss 
pY_ss pYNM_ss;

parameters kappa_ss varrho_ss varomega_ss rhox_ss z_ss a_ss yCo_ss g_ss 
Rstar_ss pistar_ss zetao_ss zetau_ss ystar_ss pCostar_ss eR_ss;

//Definitions
parameters stb_ss sg_ss sCo_ss pCoyCo_ss pYy_ss;

%--------------------------------------------------------------------------
% 2. Calibration
%--------------------------------------------------------------------------
//Calibrated parameters 
SIGMA	= 1;        // log utility (Medina and Soto, 2007)
BETTA	= 0.9995;   // quarterly subjective discount factor (steady state real interest rate of approximately 2%)
ALPHA	= 1-0.66;   // labor share of 66% (Medina and Soto, 2007)
DELTA	= 0.06/4;   // annual depreciation rate of 6% (Medina and Soto, 2007)
EPSILON_H = 11;     // steady state price markup of 10% (Medina and Soto, 2007)
EPSILON_F = 11;     // steady state price markup of 10% (Medina and Soto, 2007)
ALPHA_W	= 1;        // if 1, then wages are also indexed to A/A(-1)
O		= 0.32;     // home bias in domestic demand of 68% (imports/domestic demand=32%, 1987-2014)
CHI		= 0.61;     // CHI=c+(1-c)*t, c=0.4 (production Codelco/total, 1987-2014), t = 0.35 (general tax)
MU_CTIL	= 0;        // mean of log-normal distribution of firms' fixed costs
RHO_eR  = 0;        // persistence of the monetary policy shock

//Targeted steady state values 
stb_ss	= 0.04;        	 // average share of trade balance / GDP, 1987-2014
sg_ss	= 0.11;        	 // average share of government consumption / GDP, 1987-2014
sCo_ss	= 0.10;        	 // average share of GDP copper sector / total, 1987-2014
pi_ss	= 1.03^.25;    	 // central bank inflation target, 2001-2014
a_ss	= 1.02^.25;    	 // quarterly balanced growth path constistent with (Albagli et al., 2015)
Rstar_ss= 1.045^.25;   	 // quarterly gross foreign nominal interest rate (Fuentes and Gredig, 2008)
xi_ss	= 1.015^.25;   	 // average quarterly gross country (EMBI Chile) spread, 2001-2014
u_ss	= 0.08;        	 // average unemployment rate, 1987-2014 (=2001-2014)
e_ss	= 0.7;         	 // steady state vacancy filling probability (Cooley and Quadrini, 1999; Den Haan et al., 2000; Trigari, 2009)
rho_ss	= 0.04/(1-0.47); // Jones and Naudon (2009)
rhox_ss	= rho_ss*2/3;    // exogenous seperations 2/3 of total (Den Haan et al., 2000, and others)

//Steady state normalizations
varrho_ss	= 1;         
varomega_ss	= 1;       
z_ss		= 1;               
zetao_ss	= 1;          
zetau_ss	= 1;          
ystar_ss	= 1;          
eR_ss		= 1;             
pCostar_ss	= 1;        
pH_ss		= 1;             
h_ss		= 1;

//Pre-estimated parameters (by MLE)
RHO_g        = 0.8057;
RHO_Rstar    = 0.9677;
RHO_pistar   = 0.4488; 
RHO_ystar    = 0.8748;  
RHO_pCostar  = 0.9637;
SIG_g        = 0.0186; 
SIG_Rstar    = 0.0010; 
SIG_pistar   = 0.0259;
SIG_ystar    = 0.0054; 
SIG_pCostar  = 0.1271;

//Estimated parameters (posterior mean)
UPSILON     =	0.104341405;
PHI         =	5.094580172;
VARSIGMA    =	0.745154441;
GAMA        =	0.805587274;
ETA         =	3.525953529;
PSI         =	0.005311765;
ETASTAR     =	0.227473018;
VARPHI      =	0.268826336;
MU          =	0.747640825;
SIGMA_CTIL  =	0.30152978;
THETA_W     =	0.7787198;
VARTHETA_W  =	0.236594212;
THETA_H     =	0.328768341;
VARTHETA_H  =	0.509342846;
THETA_F     =	0.818346665;
VARTHETA_F  =	0.635874799;
RHO_R       =	0.826800084;
ALPHA_PI    =	1.545701592;
ALPHA_Y     =	0.154473528;
	
RHO_kappa	=	0.75; 0.771547576; 
RHO_varrho	=	0.75; 0.639468789; 
RHO_varomega=	0.75; 0.647259486; 
RHO_rhox	=	0.75; 0.877520566;
RHO_z       =	0.75; 0.717311389; 
RHO_a       =	0.75/2; 0.347204443; 
RHO_yCo 	=	0.75; 0.948801805; 
RHO_zetao	=	0.75; 0.873638845; 
RHO_zetau	=	0.75; 0.773949583; 
SIG_kappa	=	0.01; 0.030253667; 
SIG_varrho	=	0.01; 0.0413235; 
SIG_varomega=	0.01; 0.025738754; 
SIG_rhox	=	0.01; 0.152269856;
SIG_z       =	0.01; 0.007196123; 
SIG_a       =	0.01; 0.003037889; 
SIG_yCo	    =	0.01; 0.066067589; 
SIG_zetao	=	0.01/4; 0.01; 0.000760844;
SIG_zetau	=	0.01/4; 0.01; 0.007275765;
SIG_eR      =	0.01/4; 0.01; 0.001662937;


%--------------------------------------------------------------------------
% 3. Equilibrium equations
%--------------------------------------------------------------------------
model;
//Household (5)
exp(lam)=(exp(c)-VARSIGMA*exp(c(-1))/exp(a(-1)))^(-SIGMA); //E1
1/exp(q)=(1-GAMA/2*(exp(i)/exp(i(-1))*exp(a(-1))-a_ss)^2-GAMA*(exp(i)/exp(i(-1))*exp(a(-1))-a_ss)*exp(i)/exp(i(-1))*exp(a(-1)))*exp(varomega)+BETTA/exp(a)^SIGMA*GAMA*exp(varrho(+1))/exp(varrho)*exp(q(+1))/exp(q)*exp(lam(+1))/exp(lam)*(exp(i(+1))/exp(i)*exp(a)-a_ss)*(exp(i(+1))/exp(i)*exp(a))^2*exp(varomega(+1)); //E2
exp(q)=BETTA/exp(a)^SIGMA*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)*(exp(rK(+1))+exp(q(+1))*(1-DELTA)); //E3
exp(lam)=BETTA/exp(a)^SIGMA*exp(R)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(pi(+1)); //E4
exp(k)=(1-DELTA)*exp(k(-1))/exp(a(-1))+(1-GAMA/2*(exp(i)/exp(i(-1))*exp(a(-1))-a_ss)^2)*exp(varomega)*exp(i); //E5

//Preference shifter (2)
exp(theta)=exp(chitil)*(exp(c)-VARSIGMA*exp(c(-1))/exp(a(-1)))^-SIGMA; //E6
exp(chitil)=exp(chitil(-1))^(1-UPSILON)*(exp(c)-VARSIGMA*exp(c(-1))/exp(a(-1)))^(SIGMA*UPSILON); //E7

//Labor market (6)
exp(m) = MBAR*exp(v)^(1-MU)*u^MU; //E8
exp(n)=(1-exp(rho))*(exp(n(-1))+exp(m(-1))); //E9
u=1-exp(n); //E10
exp(s)=exp(m)/u; //E11
exp(e)=exp(m)/exp(v); //E12
exp(rho)=exp(rhox)+(1-exp(rhox))*exp(rhon); //E13
exp(thetaL)=exp(v)/u; //Labor market slackness

//Average value of employment (5)
exp(Ebar) = exp(w)*exp(h) - exp(theta)*exp(kappa)*exp(h)^(1+PHI)/(exp(lam)*(1+PHI)) 
+ exp(a)*BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)
*((1-exp(rho(+1)))*(exp(Estar(+1))+THETA_W*exp(Delbar(+1)))+exp(rho(+1))*exp(U(+1))); //E14
exp(Estar) = exp(wstar)*exp(h) - exp(theta)*exp(kappa)*exp(h)^(1+PHI)/(exp(lam)*(1+PHI)) 
+ exp(a)*BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)
*((1-exp(rho(+1)))*(exp(Estar(+1))+THETA_W*exp(Delstar(+1)))+exp(rho(+1))*exp(U(+1))); //E15
exp(U) = BBAR + exp(a)*BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)
*(exp(s)*(1-exp(rho(+1)))*exp(Ebar(+1))+(1-exp(s)*(1-exp(rho(+1))))*exp(U(+1))); //E16
exp(Delbar) = (exp(a(-1))^ALPHA_W*exp(pi(-1))^VARTHETA_W*pi_ss^(1-VARTHETA_W)/exp(a(-1))*exp(w(-1)) - exp(wstar)) * exp(h)
+ exp(a) * THETA_W * BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)*exp(Delbar(+1)); //E17
exp(Delstar) = (exp(a(-1))^ALPHA_W*exp(pi(-1))^VARTHETA_W*pi_ss^(1-VARTHETA_W)/exp(a(-1))*exp(wstar(-1)) - exp(wstar)) * exp(h)
+ exp(a) * THETA_W * BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)*exp(Delstar(+1)); //E18

//Negotiated wage (6) 
exp(epsi) = exp(h) + THETA_W * exp(a)^ALPHA_W*exp(pi)^VARTHETA_W*pi_ss^(1-VARTHETA_W) * (1-exp(rho(+1))) * BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam) 
* exp(epsi(+1))/exp(pi(+1)); //E19 
exp(mu) = exp(h) + THETA_W * exp(a)^ALPHA_W*exp(pi)^VARTHETA_W*pi_ss^(1-VARTHETA_W) * (1-exp(e))  * (1-exp(rho(+1)))* BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam) 
* exp(mu(+1))/exp(pi(+1)); //E20 
exp(Fstar) = exp(pm)*(1-ALPHA)*exp(yH)*exp(DelH)/exp(n) - exp(wstar)*exp(h) - exp(hC) + OMEGA/exp(e); //E21
exp(Wstar) = VARPHI/(1-VARPHI)* exp(epsi)/exp(mu) * exp(Fstar); //E22
exp(Wbar) = exp(Ebar) - BBAR - exp(a)*BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)*(exp(s)*(1-exp(rho(+1)))*exp(Wbar(+1))+ exp(U(+1))); //E23
exp(Wstar) = exp(wstar)*exp(h)- exp(theta)*exp(kappa)*exp(h)^(1+PHI)/(exp(lam)*(1+PHI))- BBAR+ exp(a)*BETTA/(exp(a)^SIGMA)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam) * (1-exp(rho(+1)))
*(THETA_W*VARPHI/(1-VARPHI)*exp(epsi(+1))/exp(mu(+1))*(1-exp(a(+1))^ALPHA_W*exp(pi(+1))^VARTHETA_W*pi_ss^(1-VARTHETA_W))*exp(wstar(+1))*exp(h(+1))- exp(s)*exp(Wbar(+1)) + exp(Wstar(+1))); //E24

//Average real wage (3)
exp(w) = (1-THETA_W)*exp(wstar) + THETA_W*exp(wbar); //E25
exp(wbar)*exp(pi) = exp(a(-1))^ALPHA_W*exp(pi(-1))^VARTHETA_W*pi_ss^(1-VARTHETA_W) / (exp(a(-1))*exp(n)) * (1-exp(rho))*( (1-exp(s(-1))) * exp(w(-1)) * exp(n(-1)) + exp(s(-1))*exp(wtil(-1)) ); //E26
exp(wtil)*exp(pi) = (1-THETA_W)*exp(wstar)*exp(pi) + THETA_W * (exp(a(-1))^ALPHA_W*exp(pi(-1))^VARTHETA_W*pi_ss^(1-VARTHETA_W)) / exp(a(-1)) * exp(wtil(-1)); //E27

//Individual hours
exp(h)=(exp(pm)*(1-ALPHA)^2*exp(lam)*exp(yH)*exp(DelH)/(exp(theta)*exp(kappa)*exp(n)))^(1/(1+PHI)); //E28

//Final and composite goods (5)
exp(yC)=((1-O)^(1/ETA)*exp(xH)^((ETA-1)/ETA)+O^(1/ETA)*exp(xF)^((ETA-1)/ETA))^(ETA/(ETA-1)); //E29
exp(xH)=(1-O)*exp(pH)^(-ETA)*exp(yC); //E30
exp(xF)=O*exp(pF)^(-ETA)*exp(yC); //E31
exp(yH)*exp(DelH)=exp(z)*(exp(k(-1))/exp(a(-1)))^ALPHA*(exp(a)*exp(h)*exp(n))^(1-ALPHA); //E32
exp(yF)*exp(DelF)=exp(imp); //E33

//Calvo-pricing of retail goods (9)
exp(mcH)=exp(pm)/exp(pH); //E34
exp(mcF)=exp(rer)/exp(pF); //E35
exp(fH)=exp(ptilH)^-EPSILON_H*exp(mcH)*exp(yH)+BETTA/exp(a)^(SIGMA-1)*THETA_H*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi(+1)))^-EPSILON_H*(exp(ptilH)/exp(ptilH(+1)))^-EPSILON_H*(exp(pH)/exp(pH(+1)))^(-1-EPSILON_H)*exp(fH(+1))); //E36
exp(fH)=exp(ptilH)^(1-EPSILON_H)*exp(yH)*(EPSILON_H-1)/EPSILON_H+BETTA/exp(a)^(SIGMA-1)*THETA_H*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi(+1)))^(1-EPSILON_H)*(exp(ptilH)/exp(ptilH(+1)))^(1-EPSILON_H)*(exp(pH)/exp(pH(+1)))^-EPSILON_H*exp(fH(+1))); //E37
exp(fF)=exp(ptilF)^-EPSILON_F*exp(mcF)*exp(yF)+BETTA/exp(a)^(SIGMA-1)*THETA_F*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi(+1)))^-EPSILON_F*(exp(ptilF)/exp(ptilF(+1)))^-EPSILON_F*(exp(pF)/exp(pF(+1)))^(-1-EPSILON_F)*exp(fF(+1))); //E38
exp(fF)=exp(ptilF)^(1-EPSILON_F)*exp(yF)*(EPSILON_F-1)/EPSILON_F+BETTA/exp(a)^(SIGMA-1)*THETA_F*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi(+1)))^(1-EPSILON_F)*(exp(ptilF)/exp(ptilF(+1)))^(1-EPSILON_F)*(exp(pF)/exp(pF(+1)))^-EPSILON_F*exp(fF(+1))); //E39
1=(1-THETA_H)*exp(ptilH)^(1-EPSILON_H)+THETA_H*(exp(pH(-1))/exp(pH)*exp(pi(-1))^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi))^(1-EPSILON_H); //E40
1=(1-THETA_F)*exp(ptilF)^(1-EPSILON_F)+THETA_F*(exp(pF(-1))/exp(pF)*exp(pi(-1))^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi))^(1-EPSILON_F); //E41
OMEGA/exp(e)=BETTA/exp(a)^(SIGMA-1)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)*(1-exp(rho(+1)))*(exp(pm(+1))*(1-ALPHA)*exp(yH(+1))*exp(DelH(+1))/exp(n(+1))-exp(w(+1))*exp(h(+1))-exp(hC(+1))+OMEGA/exp(e(+1))); //E61

//Intermediate goods (6)
exp(rK)=exp(pm)*ALPHA*exp(yH)*exp(a(-1))/exp(k(-1)); //E43
exp(yC)=exp(c)+exp(i)+exp(g)+exp(n)*exp(hC)+OMEGA*exp(v); //E44
exp(rhon)=1-normcdf((log(exp(ctil))-MU_CTIL)/SIGMA_CTIL); //E45
exp(ctil)=exp(pm)*(1-ALPHA)*exp(yH)*exp(DelH)/exp(n)-exp(w)*exp(h)+OMEGA/exp(e); //E46         
exp(hC)=exp(MU_CTIL+SIGMA_CTIL^2/2)*normcdf((log(exp(ctil))-MU_CTIL-SIGMA_CTIL^2)/SIGMA_CTIL)/(1-exp(rhon)); //E47
exp(y)=exp(c)+exp(i)+exp(g)+exp(xHstar)+exp(yCo)-exp(imp); //E48

//Taylor rule (1)
exp(R)/R_ss=(exp(R(-1))/R_ss)^RHO_R*((exp(pi)/pi_ss)^ALPHA_PI*(exp(y)/exp(y(-1)))^ALPHA_Y)^(1-RHO_R)*exp(eR); //E49

//Rest of the world (4)
exp(xi)=xi_ss*exp(PSI*(exp(rer)*dstar-rer_ss*dstar_ss)/(rer_ss*dstar_ss)+(exp(zetao)-zetao_ss)/zetao_ss+(exp(zetau)-zetau_ss)/zetau_ss); //E50
exp(rer)/exp(rer(-1))=exp(piS)*exp(pistar)/exp(pi); //E51
exp(lam)=BETTA/exp(a)^SIGMA*exp(Rstar)*exp(xi)*exp(varrho(+1))/exp(varrho)*exp(piS(+1))*exp(lam(+1))/exp(pi(+1)); //E52
exp(xHstar)=OSTAR*(exp(pH)/exp(rer))^-ETASTAR*exp(ystar); //E53

//Aggregation and market clearing (7)
exp(yH)=exp(xH)+exp(xHstar); //E54
exp(yF)=exp(xF); //E55
exp(DelH)=(1-THETA_H)*exp(ptilH)^-EPSILON_H+THETA_H*(exp(pH(-1))/exp(pH)*exp(pi(-1))^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi))^-EPSILON_H*exp(DelH(-1)); //E56
exp(DelF)=(1-THETA_F)*exp(ptilF)^-EPSILON_F+THETA_F*(exp(pF(-1))/exp(pF)*exp(pi(-1))^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi))^-EPSILON_F*exp(DelF(-1)); //E57
exp(rer)*dstar=exp(rer)*dstar(-1)/exp(a(-1))/exp(pistar)*exp(Rstar(-1))*exp(xi(-1))-tb+(1-CHI)*exp(rer)*exp(pCostar)*exp(yCo); //E58
tb=exp(pH)*exp(xHstar)-exp(rer)*exp(imp)+exp(rer)*exp(pCostar)*exp(yCo); //E59
exp(pY)*exp(y)=exp(c)+exp(i)+exp(g)+tb; //E60

// Exogenous AR(1) processes (15)
kappa-log(kappa_ss)=RHO_kappa*(kappa(-1)-log(kappa_ss))+eps_kappa; 
varrho-log(varrho_ss)=RHO_varrho*(varrho(-1)-log(varrho_ss))+eps_varrho; 
varomega-log(varomega_ss)=RHO_varomega*(varomega(-1)-log(varomega_ss))+eps_varomega; 
rhox-log(rhox_ss)=RHO_rhox*(rhox(-1)-log(rhox_ss))+eps_rhox; 
z-log(z_ss)=RHO_z*(z(-1)-log(z_ss))+eps_z; 
a-log(a_ss)=RHO_a*(a(-1)-log(a_ss))+eps_a; 
yCo-log(yCo_ss)=RHO_yCo*(yCo(-1)-log(yCo_ss))+eps_yCo; 
g-log(g_ss)=RHO_g*(g(-1)-log(g_ss))+eps_g; 
zetao-log(zetao_ss)=RHO_zetao*(zetao(-1)-log(zetao_ss))+eps_zetao; 
zetau-log(zetau_ss)=RHO_zetau*(zetau(-1)-log(zetau_ss))+eps_zetau; 
Rstar-log(Rstar_ss)=RHO_Rstar*(Rstar(-1)-log(Rstar_ss))+eps_Rstar; 
pistar-log(pistar_ss)=RHO_pistar*(pistar(-1)-log(pistar_ss))+eps_pistar;
ystar-log(ystar_ss)=RHO_ystar*(ystar(-1)-log(ystar_ss))+eps_ystar;
pCostar-log(pCostar_ss)=RHO_pCostar*(pCostar(-1)-log(pCostar_ss))+eps_pCostar; 
eR-log(eR_ss)=RHO_eR*(eR(-1)-log(eR_ss))+eps_eR; 

//Definitions (4)
exp(y) = exp(yNM) + exp(yCo);
exp(pYNM)*exp(yNM) = exp(pY)*exp(y)- exp(pF)*exp(pCostar)*exp(yCo);
stb=tb/(exp(pY)*exp(y));
exp(gam_Y)=exp(pY)*exp(y)/(exp(pY(-1))*exp(y(-1)))*exp(a(-1));
exp(gam_YNM)=exp(pYNM)*exp(yNM)/(exp(pYNM(-1))*exp(yNM(-1)))*exp(a(-1));
exp(gam_C)=exp(c)/exp(c(-1))*exp(a(-1)); 
exp(gam_I)=exp(i)/exp(i(-1))*exp(a(-1)); 
exp(gam_W)=exp(w)/exp(w(-1))*exp(a(-1));
exp(pCoyCo) = exp(rer)*exp(pCostar)*exp(yCo);

//Measurement equations (16)
gam_Y_obs=100*(gam_Y-log(a_ss)); 
gam_YNM_obs=100*(gam_YNM-log(a_ss)); 
gam_C_obs=100*(gam_C-log(a_ss)); 
gam_I_obs=100*(gam_I-log(a_ss)); 
gam_W_obs=100*(gam_W-log(a_ss)); 
R_obs=100*(R-log(R_ss)); 
xi_obs=100*(xi-log(xi_ss)-(zetau-log(zetau_ss))); 
pi_obs=100*(pi-log(pi_ss)); 
rer_obs=100*(rer-log(rer_ss)); 
hn_obs=100*(h-log(h_ss)+n-log(n_ss)); 
u_obs=100*(u-u_ss); 
pCoyCo_obs=100*(pCoyCo-log(pCoyCo_ss));
g_obs=100*(g-log(g_ss)); 
Rstar_obs=100*(Rstar-log(Rstar_ss)); 
pistar_obs=100*(pistar-log(pistar_ss)); 
ystar_obs=100*(ystar-log(ystar_ss)); 
pCostar_obs=100*(pCostar-log(pCostar_ss)); 

end;

%--------------------------------------------------------------------------
% 4. Steady state
%--------------------------------------------------------------------------
steady_state_model;

//Computing the steady state and endogenous parameters
theta_ss=1;
R_ss=a_ss^SIGMA*pi_ss/BETTA;
q_ss=1/varomega_ss;
rK_ss=q_ss*(a_ss^SIGMA/BETTA-1+DELTA);
piS_ss=a_ss^SIGMA*pi_ss/(BETTA*Rstar_ss*xi_ss);
pistar_ss=pi_ss/piS_ss;
ptilH_ss=1;
ptilF_ss=1;
DelH_ss=ptilH_ss^-EPSILON_H;
DelF_ss=ptilF_ss^-EPSILON_F;
mcH_ss=(EPSILON_H-1)/EPSILON_H*ptilH_ss;
mcF_ss=(EPSILON_F-1)/EPSILON_F*ptilF_ss;
pm_ss=pH_ss*mcH_ss;
n_ss=1-u_ss;
rhon_ss=(rho_ss-rhox_ss)/(1-rhox_ss);
ctil_ss=exp(MU_CTIL+SIGMA_CTIL*norminv(1-rhon_ss));        
hC_ss=exp(MU_CTIL+SIGMA_CTIL^2/2)*normcdf((log(ctil_ss)-MU_CTIL-SIGMA_CTIL^2)/SIGMA_CTIL)/(1-rhon_ss);
OMEGA=e_ss*a_ss^(1-SIGMA)*BETTA*(1-rho_ss)*(ctil_ss-hC_ss);
v_ss=rho_ss*n_ss/(e_ss*(1-rho_ss));
m_ss=e_ss*v_ss;
MBAR=m_ss*v_ss^(MU-1)*u_ss^-MU; 
s_ss=m_ss/u_ss; 		
k_ss=a_ss^2*h_ss*n_ss*(ALPHA*pm_ss*z_ss/(DelH_ss*rK_ss))^(1/(1-ALPHA));
yH_ss=z_ss*(k_ss/a_ss)^ALPHA*(a_ss*h_ss*n_ss)^(1-ALPHA)/DelH_ss;
thetaL_ss=v_ss/u_ss; //Labor market slackness
w_ss=(1/h_ss)*(pm_ss*(1-ALPHA)*yH_ss*DelH_ss/n_ss+OMEGA/e_ss-ctil_ss);
wstar_ss=w_ss*(1-(1-rho_ss)*(1-s_ss)*a_ss^(ALPHA_W-1)*THETA_W)/(1-THETA_W)*(1+(1-rho_ss)*s_ss*a_ss^(ALPHA_W-1)*THETA_W/(n_ss*(1-a_ss^(ALPHA_W-1)*THETA_W)))^-1;
wtil_ss=(1-THETA_W)/(1-a_ss^(ALPHA_W-1)*THETA_W)*wstar_ss;
wbar_ss=w_ss/THETA_W-(1-THETA_W)/THETA_W*wstar_ss;
epsi_ss=h_ss/(1-THETA_W*a_ss^(ALPHA_W-SIGMA)*BETTA*(1-rho_ss));
mu_ss=h_ss/(1-THETA_W*a_ss^(ALPHA_W-SIGMA)*BETTA*(1-rho_ss)*(1-e_ss));
Fstar_ss=pm_ss*(1-ALPHA)*yH_ss*DelH_ss/n_ss-wstar_ss*h_ss-hC_ss+OMEGA/e_ss;
Wstar_ss=VARPHI/(1-VARPHI)*epsi_ss/mu_ss*Fstar_ss;
Delbar_ss =(a_ss^(ALPHA_W-1)*pi_ss*w_ss/wstar_ss-1)/(1-a_ss^(1-SIGMA)*BETTA*THETA_W)*wstar_ss*h_ss;
Delstar_ss=(a_ss^(ALPHA_W-1)*pi_ss-1)/(1-a_ss^(1-SIGMA)*BETTA*THETA_W)*wstar_ss*h_ss;
AUX1=pm_ss*(1-ALPHA)^2*yH_ss*DelH_ss/(theta_ss*h_ss^(1+PHI)*n_ss);
alpha1=(wstar_ss*h_ss- AUX1*theta_ss*h_ss^(1+PHI)/(1+PHI)+ a_ss^(1-SIGMA)*BETTA*(1-rho_ss)*THETA_W*Delstar_ss)/(1-a_ss^(1-SIGMA)*BETTA*(1-rho_ss));
alpha2=a_ss^(1-SIGMA)*BETTA*rho_ss/(1-a_ss^(1-SIGMA)*BETTA*(1-rho_ss));
alpha3=w_ss*h_ss-AUX1*theta_ss*h_ss^(1+PHI)/(1+PHI)+a_ss^(1-SIGMA)*BETTA*(1-rho_ss)*(THETA_W*Delbar_ss+alpha1);
alpha4=((1-rho_ss)*alpha2+rho_ss)*a_ss^(1-SIGMA)*BETTA;
alpha5=-a_ss^(1-SIGMA)*BETTA*s_ss*(1-rho_ss)*alpha3;
alpha6=1-a_ss^(1-SIGMA)*BETTA*(1+s_ss*(1-rho_ss)*(alpha4-1));
alpha7=((1+a_ss^(1-SIGMA)*BETTA*(1-rho_ss)*THETA_W*VARPHI*epsi_ss/((1-VARPHI)*mu_ss)*(1-a_ss^ALPHA_W*pi_ss))*wstar_ss*h_ss- AUX1*theta_ss*h_ss^(1+PHI)/(1+PHI) - (1-a_ss^(1-SIGMA)*BETTA*(1-rho_ss))*Wstar_ss - alpha5)/(a_ss^(1-SIGMA)*BETTA*s_ss*(1-rho_ss));
alpha8=alpha6/(a_ss^(1-SIGMA)*BETTA*s_ss*(1-rho_ss));
U_ss=(alpha7-alpha3+alpha5+a_ss^(1-SIGMA)*BETTA*s_ss*(1-rho_ss)*alpha7)/(alpha8+alpha4-alpha6+a_ss^(1-SIGMA)*BETTA*s_ss*(1-rho_ss)*alpha8-a_ss^(1-SIGMA)*BETTA);
Estar_ss=alpha1+alpha2*U_ss;
Ebar_ss=alpha3+alpha4*U_ss;
BBAR=alpha5+alpha6*U_ss;
Wbar_ss=alpha7-alpha8*U_ss;
fH_ss=ptilH_ss^-EPSILON_H*yH_ss*mcH_ss/(1-BETTA*a_ss^(1-SIGMA)*THETA_H);
i_ss=k_ss*(1-(1-DELTA)/a_ss)/varomega_ss;
pF_ss=(1/O-(1-O)/O*pH_ss^(1-ETA))^(1/(1-ETA));
rer_ss=mcF_ss*pF_ss;
pYy_ss=(pH_ss*yH_ss+(hC_ss*n_ss+OMEGA*v_ss)*(pF_ss*O*pF_ss^-ETA*(1-mcF_ss*DelF_ss)-1))/(1-sCo_ss-pF_ss*(1-mcF_ss*DelF_ss)*O*pF_ss^-ETA*(1-stb_ss));
tb_ss=stb_ss*pYy_ss;
g_ss=sg_ss*pYy_ss;
yCo_ss=sCo_ss*pYy_ss/(rer_ss*pCostar_ss);
pCoyCo_ss=rer_ss*pCostar_ss*yCo_ss;
xHstar_ss=yH_ss-(1-O)*pH_ss^-ETA*(pYy_ss-tb_ss+n_ss*hC_ss+OMEGA*v_ss);
xF_ss=(pH_ss*xHstar_ss+rer_ss*pCostar_ss*yCo_ss-tb_ss)/rer_ss;
yC_ss=(xF_ss/O)*pF_ss^ETA;
c_ss=yC_ss-g_ss-i_ss-n_ss*hC_ss-OMEGA*v_ss;
lam_ss=(c_ss-VARSIGMA*c_ss/a_ss)^(-SIGMA);
kappa_ss=pm_ss*lam_ss*(1-ALPHA)^2*yH_ss/(theta_ss*h_ss^(1+PHI)*n_ss);
chitil_ss=c_ss^SIGMA*(1-VARSIGMA/a_ss)^SIGMA;
OSTAR=xHstar_ss/ystar_ss*(pH_ss/rer_ss)^ETASTAR;
dstar_ss=-(tb_ss-(1-CHI)*rer_ss*pCostar_ss*yCo_ss)/(rer_ss*(1-Rstar_ss*xi_ss/pistar_ss/a_ss));
xH_ss=yH_ss-xHstar_ss;
yF_ss=xF_ss;
fF_ss=ptilF_ss^-EPSILON_F*yF_ss*mcF_ss/(1-BETTA*a_ss^(1-SIGMA)*THETA_F);
imp_ss=yF_ss*DelF_ss;
y_ss=c_ss+i_ss+g_ss+xHstar_ss+yCo_ss-imp_ss;
pY_ss=(c_ss+i_ss+g_ss+tb_ss)/y_ss;
yNM_ss = y_ss -yCo_ss;
pYNM_ss = (pY_ss*y_ss - pF_ss*pCostar_ss*yCo_ss)/yNM_ss;

//Initial values for Dynare numerical solver
lam=log(lam_ss);
q=log(q_ss); 
c=log(c_ss);
i=log(i_ss);
k=log(k_ss);
theta=log(theta_ss); 
chitil=log(chitil_ss);
u=(u_ss); 
thetaL=log(thetaL_ss); // Labor market slackness
m=log(m_ss); 
n=log(n_ss); 
rho=log(rho_ss); 
s=log(s_ss); 
e=log(e_ss);
Ebar=log(Ebar_ss); 
Estar=log(Estar_ss); 
U=log(U_ss); 
Delbar=log(Delbar_ss); 
Delstar=log(Delstar_ss);
wstar=log(wstar_ss); 
epsi=log(epsi_ss); 
mu=log(mu_ss); 
Fstar=log(Fstar_ss); 
Wstar=log(Wstar_ss); 
Wbar=log(Wbar_ss);
w=log(w_ss); 
wbar=log(wbar_ss); 
wtil=log(wtil_ss);
h=log(h_ss);
yC=log(yC_ss); 
xH=log(xH_ss); 
xF=log(xF_ss); 
yH=log(yH_ss); 
yF=log(yF_ss);
pm=log(pm_ss); 
ptilH=log(ptilH_ss); 
ptilF=log(ptilF_ss); 
fH=log(fH_ss); 
fF=log(fF_ss);  
mcH=log(mcH_ss);  
mcF=log(mcF_ss); 
rK=log(rK_ss);  
v=log(v_ss);  
rhon=log(rhon_ss);  
ctil=log(ctil_ss);  
hC=log(hC_ss);  
imp=log(imp_ss); 
pi=log(pi_ss); 
R=log(R_ss); 
xi=log(xi_ss);  
rer=log(rer_ss);  
piS=log(piS_ss);  
xHstar=log(xHstar_ss); 
pH=log(pH_ss);  
pF=log(pF_ss);  
DelH=log(DelH_ss);  
DelF=log(DelF_ss);  
dstar=(dstar_ss);  
tb=(tb_ss);  
y=log(y_ss);  
pY=log(pY_ss);
yNM=log(yNM_ss);
pYNM=log(pYNM_ss);
kappa=log(kappa_ss);  
varrho=log(varrho_ss);  
varomega=log(varomega_ss);  
rhox=log(rhox_ss); 
z=log(z_ss);  
a=log(a_ss);  
yCo=log(yCo_ss);  
g=log(g_ss);  
Rstar=log(Rstar_ss);  
pistar=log(pistar_ss);  
zetao=log(zetao_ss);  
zetau=log(zetau_ss);  
ystar=log(ystar_ss);  
pCostar=log(pCostar_ss);  
eR=log(eR_ss); 
stb=stb_ss;
gam_Y=log(a_ss);
gam_YNM=log(a_ss);
gam_C=log(a_ss);
gam_I=log(a_ss);
gam_W=log(a_ss); 
pCoyCo=log(pCoyCo_ss);
gam_Y_obs=0;
gam_YNM_obs=0;
gam_C_obs=0;
gam_I_obs=0;
gam_W_obs=0;
R_obs=0; 
xi_obs=0; 
pi_obs=0; 
rer_obs=0; 
hn_obs=0; 
u_obs=0;
pCoyCo_obs=0; 
g_obs=0; 
Rstar_obs=0; 
pistar_obs=0; 
ystar_obs=0; 
pCostar_obs=0;

end;

//Call steady state solver
steady;

//Check Blanchard-Kahn conditions
check;

%--------------------------------------------------------------------------
% 5. Estimation and simulation
%--------------------------------------------------------------------------
//Observed variables
varobs gam_YNM_obs gam_C_obs gam_I_obs gam_W_obs R_obs xi_obs pi_obs rer_obs
u_obs hn_obs pCoyCo_obs g_obs Rstar_obs pistar_obs ystar_obs pCostar_obs;

//Structural shocks
shocks;
var eps_kappa;    stderr SIG_kappa;
var eps_varrho;   stderr SIG_varrho;
var eps_varomega; stderr SIG_varomega;
var eps_rhox;     stderr SIG_rhox;
var eps_z;        stderr SIG_z;
var eps_a; 		  stderr SIG_a;
var eps_yCo;      stderr SIG_yCo;
var eps_g;        stderr SIG_g;
var eps_Rstar;    stderr SIG_Rstar;
var eps_pistar;   stderr SIG_pistar;
var eps_zetao;    stderr SIG_zetao;
var eps_zetau;    stderr SIG_zetau;
var eps_ystar;    stderr SIG_ystar;
var eps_pCostar;  stderr SIG_pCostar;
var eps_eR;       stderr SIG_eR;
var eps_aux;      stderr eps;

//Measurement errors (10% of the observed variance)
var gam_YNM_obs; 	stderr sqrt(1.1437^2*0.1);   
var gam_C_obs; 	stderr sqrt(1.2061^2*0.1);
var gam_I_obs; 	stderr sqrt(3.7542^2*0.1);
var gam_W_obs; 	stderr sqrt(0.5847^2*0.1);      
var R_obs; 		stderr sqrt(0.4164^2*0.0);
var xi_obs; 	stderr sqrt(0.1467^2*0.1);
var pi_obs; 	stderr sqrt(0.6916^2*0.1);
var rer_obs; 	stderr sqrt(5.1679^2*0.1);  
var u_obs; 		stderr sqrt(1.4293^2*0.1); 
var hn_obs; 	stderr sqrt(1.8690^2*0.1);
var pCoyCo_obs;     stderr sqrt(41.34765248^2*0.1);
end;

//Parameters to be estimated (already done)
estimated_params;
//Deep parameters     
UPSILON,    beta_pdf, 0.5, 0.25;
PHI, , 0, , normal_pdf, 2, 2;
VARSIGMA,   beta_pdf, 0.7, 0.1;
GAMA,       normal_pdf, 4, 1.5;
ETA,        inv_gamma_pdf, 1.5, 0.25;
PSI,        inv_gamma_pdf, 0.01, inf;
ETASTAR,    inv_gamma_pdf, 0.25, 0.1;
VARPHI,   beta_pdf, 0.5, 0.1;
MU,         beta_pdf, 0.5, 0.1;
SIGMA_CTIL, , 0, , normal_pdf, 0.1, 0.05;
THETA_W,    beta_pdf, 0.5, 0.1; //0.75, 0.1;
VARTHETA_W, beta_pdf, 0.5, 0.15;
THETA_H,    beta_pdf, 0.5, 0.1; //0.75, 0.1;
VARTHETA_H, beta_pdf, 0.5, 0.15;
THETA_F,    beta_pdf, 0.5, 0.1; //0.75, 0.1;
VARTHETA_F, beta_pdf, 0.5, 0.15;
RHO_R,      beta_pdf, 0.75, 0.1;
ALPHA_PI, , 1, , normal_pdf, 1.5, 0.1;
ALPHA_Y, , 0, ,  normal_pdf, 0.5/4, 0.05;

//Persistence of the AR(1) process
RHO_kappa,    beta_pdf, 0.75, 0.1;
RHO_varrho,   beta_pdf, 0.75, 0.1;
RHO_varomega, beta_pdf, 0.75, 0.1;
RHO_rhox,     beta_pdf, 0.75, 0.1;
RHO_z,        beta_pdf, 0.75, 0.1;
RHO_a,        beta_pdf, 0.75/2, 0.1;
RHO_yCo,      beta_pdf, 0.75, 0.1;
RHO_zetao,    beta_pdf, 0.75, 0.1;
RHO_zetau,    beta_pdf, 0.75, 0.1;

//Std. deviation of the error term of the AR(1) process
stderr eps_kappa,    inv_gamma_pdf, 0.01, inf;
stderr eps_varrho,   inv_gamma_pdf, 0.01, inf;
stderr eps_varomega, inv_gamma_pdf, 0.01, inf;
stderr eps_rhox,     inv_gamma_pdf, 0.01, inf;
stderr eps_z,        inv_gamma_pdf, 0.01, inf;
stderr eps_a,        inv_gamma_pdf, 0.01, inf;
stderr eps_yCo,      inv_gamma_pdf, 0.01, inf;
stderr eps_zetao,    inv_gamma_pdf, 0.01/4, inf;
stderr eps_zetau,    inv_gamma_pdf, 0.01/4, inf;
stderr eps_eR,       inv_gamma_pdf, 0.01/4, inf;

end;

//The following computes the posterior mode
//estimation(datafile=estim_data, xls_sheet=estim_data, mode_compute=4, mh_replic= 0 plot_priors=1); 

//The following runs MCMC replications to estimate posterior distribution. It uses the posterior mode as initial point
//estimation(datafile=estim_data, xls_sheet=estim_data, mode_compute=0, mode_file=search_mode, mh_replic=0, mh_jscale=0.35, mh_nblocks=1, plot_priors=0);

//The following simulates the model by first-order perturbation
stoch_simul(order=1,periods=0,irf=40,nograph);

//The following computes recursive forecasts
/*
@#for n in 22:59
    estimation(datafile=estim_data,xls_sheet=estim_data,plot_priors=0,first_obs=1,nobs=@{n}, mh_replic=0,forecast=10,nograph) gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs;
    save Standar_@{n}.mat       
@#endfor
*/