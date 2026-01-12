% Standard NKSOE model with sticky wages a la Schmitt-Grohé and Uribe (2006)

%--------------------------------------------------------------------------
% 0. Housekeeping (close all graphic windows) and options
%--------------------------------------------------------------------------
close all;

%--------------------------------------------------------------------------
% 1. Preamble
%--------------------------------------------------------------------------
//Endogenous variables
var pCoyCo lam c h w i k rK q y yNM yC yF yH xF xH xHstar R xi pi rer pH 
ptilH pF ptilF pY pYNM piS mcH fH DelH mcF fF DelF dstar imp tb stb sdstar 
gam_Y gam_YNM gam_C gam_I gam_W gam_K gam_YC gam_YF gam_YH gam_XF gam_XH 
gam_XHstar gam_M gam_LAM piH piF gammaW Theta chitil hd wtil mcW fW DelW; 

// Exogenous state variables
var varrho kappa varomega z a zetao zetau eR yCo Rstar pistar pCostar ystar g;

// Observed variables
var pCoyCo_obs yCo_obs gam_YNM_obs gam_Y_obs gam_C_obs gam_I_obs gam_W_obs 
R_obs xi_obs pi_obs rer_obs g_obs Rstar_obs ystar_obs pistar_obs pCostar_obs 
stb_obs hn_obs pi4_obs piw4_obs piw_obs; 

// Exogenous innovations
varexo eps_varrho eps_kappa eps_varomega eps_z eps_a eps_zetao eps_zetau 
eps_eR eps_yCo eps_Rstar eps_pistar eps_pCostar eps_ystar eps_g eps_aux 
eps_aux1;

// Parameters    
parameters SIGMA PHI ALPHA DELTA EPSILON_H EPSILON_F EPSILON_W O CHI UPSILON;
parameters VARSIGMA PSI ETA ETASTAR VARKAPPA_X GAMA THETA_W VARTHETA_W 
ALPHA_W THETA_H VARTHETA_H THETA_F VARTHETA_F RHO_R ALPHA_PI ALPHA_Y;
parameters BETTA OSTAR;
parameters RHO_varrho RHO_kappa RHO_varomega RHO_z RHO_a RHO_zetao RHO_zetau 
RHO_eR RHO_yCo RHO_Rstar RHO_pistar RHO_pCostar RHO_ystar RHO_g;
parameters SIG_varrho SIG_kappa SIG_varomega SIG_z SIG_a SIG_zetao SIG_zetau 
SIG_eR SIG_yCo SIG_Rstar SIG_pistar SIG_pCostar SIG_ystar SIG_g;
parameters pCoyCo_ss lam_ss c_ss h_ss w_ss i_ss k_ss rK_ss q_ss y_ss yNM_ss 
yC_ss yF_ss yH_ss xF_ss xH_ss xHstar_ss R_ss xi_ss pi_ss rer_ss piS_ss pH_ss 
ptilH_ss pF_ss ptilF_ss pY_ss pm_ss mcH_ss fH_ss DelH_ss mcF_ss fF_ss DelF_ss 
dstar_ss imp_ss tb_ss stb_ss sdstar_ss gammaW_ss Theta_ss chitil_ss;
parameters hd_ss wtil_ss mcW_ss fW_ss DelW_ss;
parameters varrho_ss kappa_ss varomega_ss z_ss a_ss zetao_ss zetau_ss eR_ss 
yCo_ss Rstar_ss pistar_ss pCostar_ss ystar_ss g_ss;
parameters sCo_ss sg_ss pYy_ss;

//Calibrated parameters (Endogenous BETTA, kappa_ss, OMEGA, pistar_ss, ALPHA, g_ss, yCo_ss)
SIGMA=1;         // log utility (Medina and Soto, 2007)
PHI=1;               // unitary Frisch elasticity (Adolfson et al., 2008)
ALPHA=1-0.66;        // labor share of 66% (Medina and Soto, 2007)
DELTA=0.06/4;        // annual depreciation rate of 6% (Medina and Soto, 2007)
EPSILON_H=11;        // steady state price markup of 10% (Medina and Soto, 2007)
EPSILON_F=11;        // steady state price markup of 10% (Medina and Soto, 2007)
EPSILON_W=11;        // steady state wage markup of 10% (Medina and Soto, 2007)
ALPHA_W=1;       // 0 for no indexation to A/A(-1), 1 for indexation
O=0.32;              // home bias in domestic demand of 68% (imports/domestic demand=32%, 1987-2014)
CHI=0.61;            // CHI=c+(1-c)*t, c=0.4 (production Codelco/total, 1987-2014), t = 0.35 (general tax)
VARKAPPA_X=0;        
RHO_eR=0;            // Monetary policy persistence

// Targeted steady state values (watch out, sdstar_ss should be around 2 to be constistent with NFA/GDP ratio of 50%) 
stb_ss=0.04;         // average share of trade balance / GDP, 1987-2014
sg_ss=0.11;          // average share of government consumption / GDP, 1987-2014
sCo_ss=0.10;         // average share of GDP copper sector / total, 1987-2014
pi_ss=1.03^.25;      // central bank inflation target, 2001-2014
a_ss=1.02^.25;       // quarterly balanced growth path constistent with (Albagli et al., 2015)
BETTA=0.9995;        // quarterly subjective discount factor (steady state real interest rate of approximately 2%)
Rstar_ss=1.045^.25;  // quarterly gross foreign nominal interest rate (Fuentes and Gredig, 2008)
xi_ss=1.015^.25;     // average quarterly gross country (EMBI Chile) spread, 2001-2014
   
// steady state normalizations
h_ss=0.3;            
varrho_ss=1;         
varomega_ss=1;       
pH_ss=1;             
z_ss=1;               
zetao_ss=1;          
zetau_ss=1;          
eR_ss=1;             
ystar_ss=1;          
pCostar_ss=1;        

// Pre-estimated parameters (exogenous processes)
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

// Estimated parameters (posterior mean)
UPSILON      =   0.341650151;
PHI          =   1.002834253; 
VARSIGMA     =   0.817695941; 
PSI          =   0.004465966; 
ETA          =   1.736828268; 
ETASTAR      =   0.180768964; 
GAMA         =   5.479084145; 
THETA_W      =   0.900460184; 
VARTHETA_W   =   0.402842622; 
THETA_H      =   0.555152329; 
VARTHETA_H   =   0.213842068; 
THETA_F      =   0.56949997; 
VARTHETA_F   =   0.542133562; 
RHO_R        =   0.853519042; 
ALPHA_PI     =   1.498534036; 
ALPHA_Y      =   0.128096346; 

RHO_yCo         =	0.75; 0.94410615; 
RHO_varrho      =	0.75; 0.73967766; 
RHO_kappa       =	0.75; 0.736565274; 
RHO_varomega	=	0.75; 0.845410671; 
RHO_z           =	0.75; 0.497756816; 
RHO_a           =	0.75/2; 0.316593822; 
RHO_zetao       =	0.75; 0.849218072; 
RHO_zetau       =	0.75; 0.743732377; 
SIG_yCo         =	0.01; 0.06847054; 
SIG_varrho      =	0.01; 0.06121832; 
SIG_kappa       =	0.01; 0.127529273; 
SIG_varomega	=	0.01; 0.074980488; 
SIG_z           =	0.01; 0.0186898; 
SIG_a           =	0.01; 0.018796754; 
SIG_zetao       =	0.01/4; 0.00073068; 
SIG_zetau       =	0.01/4; 0.008726603; 
SIG_eR          =	0.01/4; 0.001612481; 

%--------------------------------------------------------------------------
% 2. Model (89=41+17+16+15 equations)
%--------------------------------------------------------------------------
model; 
// Equilibrium conditions 
exp(lam)=(exp(c)-VARSIGMA*exp(c(-1))/exp(a(-1)))^(-SIGMA); //E1
exp(w)*exp(mcW)=exp(Theta)*exp(kappa)*exp(h)^PHI/exp(lam); //E2
exp(lam)=BETTA/exp(a)^SIGMA*exp(R)*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(pi(+1)); //E3
exp(lam)=BETTA/exp(a)^SIGMA*exp(Rstar)*exp(xi)*exp(varrho(+1))/exp(varrho)*exp(piS(+1))*exp(lam(+1))/exp(pi(+1)); //E4
exp(yC)=((1-O)^(1/ETA)*exp(xH)^((ETA-1)/ETA)+O^(1/ETA)*exp(xF)^((ETA-1)/ETA))^(ETA/(ETA-1)); //E5
exp(xF)=O*exp(pF)^(-ETA)*exp(yC); //E6
exp(xH)=(1-O)*exp(pH)^(-ETA)*exp(yC); //E7
exp(mcH)=exp(rK)^ALPHA*exp(w)^(1-ALPHA)/(exp(pH)*exp(z)*exp(a)^(1-ALPHA))/(ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA)); //E8
exp(fH)=exp(ptilH)^-EPSILON_H*exp(yH)*exp(mcH)+BETTA/exp(a)^(SIGMA-1)*THETA_H*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi(+1)))^-EPSILON_H*(exp(ptilH)/exp(ptilH(+1)))^-EPSILON_H*(exp(pH)/exp(pH(+1)))^(-1-EPSILON_H)*exp(fH(+1))); //E9
exp(fH)=exp(ptilH)^(1-EPSILON_H)*exp(yH)*(EPSILON_H-1)/EPSILON_H+BETTA/exp(a)^(SIGMA-1)*THETA_H*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi(+1)))^(1-EPSILON_H)*(exp(ptilH)/exp(ptilH(+1)))^(1-EPSILON_H)*(exp(pH)/exp(pH(+1)))^-EPSILON_H*exp(fH(+1))); //E10
exp(xHstar)=VARKAPPA_X*exp(xHstar(-1))/exp(a(-1))+(1-VARKAPPA_X)*OSTAR*(exp(pH)/exp(rer))^-ETASTAR*exp(ystar); //E11
exp(R)/R_ss=(exp(R(-1))/R_ss)^RHO_R*((exp(pi)/pi_ss)^ALPHA_PI*(exp(y)/exp(y(-1)))^ALPHA_Y)^(1-RHO_R)*exp(eR); //E12
exp(yH)*exp(DelH)=exp(z)*(exp(k(-1))/exp(a(-1)))^ALPHA*(exp(a)*exp(hd))^(1-ALPHA); //E13
1=(1-THETA_H)*exp(ptilH)^(1-EPSILON_H)+THETA_H*(exp(pH(-1))/exp(pH)*exp(pi(-1))^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi))^(1-EPSILON_H); //E14
exp(yH)=exp(xH)+exp(xHstar); //E15
exp(yC)=exp(c)+exp(i)+exp(g); //E16
exp(rer)/exp(rer(-1))=exp(piS)*exp(pistar)/exp(pi); //E17
exp(y)=exp(c)+exp(i)+exp(g)+exp(xHstar)+exp(yCo)-exp(imp); //E18
tb=exp(pH)*exp(xHstar)-exp(rer)*exp(imp)+exp(rer)*exp(pCostar)*exp(yCo); //E19
exp(rer)*dstar=exp(rer)*dstar(-1)/exp(a(-1))/exp(pistar)*exp(Rstar(-1))*exp(xi(-1))-tb+(1-CHI)*exp(rer)*exp(pCostar)*exp(yCo); //E20
exp(xi)=xi_ss*exp(PSI*(exp(rer)*dstar-rer_ss*dstar_ss)/(rer_ss*dstar_ss)+(exp(zetao)-zetao_ss)/zetao_ss+(exp(zetau)-zetau_ss)/zetau_ss); //E21, for calibrations with dstar_ss<0, write exp(PSI*(exp(rer)*dstar-rer_ss*dstar_ss)/(-rer_ss*dstar_ss)...
exp(DelH)=(1-THETA_H)*exp(ptilH)^-EPSILON_H+THETA_H*(exp(pH(-1))/exp(pH)*exp(pi(-1))^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pi))^-EPSILON_H*exp(DelH(-1)); //E22
exp(k)=(1-DELTA)*exp(k(-1))/exp(a(-1))+(1-GAMA/2*(exp(i)/exp(i(-1))*exp(a(-1))-a_ss)^2)*exp(varomega)*exp(i); //E23
exp(q)=BETTA/exp(a)^SIGMA*exp(varrho(+1))/exp(varrho)*exp(lam(+1))/exp(lam)*(exp(rK(+1))+exp(q(+1))*(1-DELTA)); //E24
exp(k(-1))/exp(hd)=exp(a(-1))*ALPHA/(1-ALPHA)*exp(w)/exp(rK); //E25
1/exp(q)=(1-GAMA/2*(exp(i)/exp(i(-1))*exp(a(-1))-a_ss)^2-GAMA*(exp(i)/exp(i(-1))*exp(a(-1))-a_ss)*exp(i)/exp(i(-1))*exp(a(-1)))*exp(varomega)+BETTA/exp(a)^SIGMA*GAMA*exp(varrho(+1))/exp(varrho)*exp(q(+1))/exp(q)*exp(lam(+1))/exp(lam)*(exp(i(+1))/exp(i)*exp(a)-a_ss)*(exp(i(+1))/exp(i)*exp(a))^2*exp(varomega(+1)); //E26
exp(pY)*exp(y)=exp(c)+exp(i)+exp(g)+tb; //E27
1=(1-THETA_F)*exp(ptilF)^(1-EPSILON_F)+THETA_F*(exp(pF(-1))/exp(pF)*exp(pi(-1))^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi))^(1-EPSILON_F); //E28
exp(yF)=exp(xF); //E29
exp(imp)=exp(yF)*exp(DelF); //E30
exp(DelF)=(1-THETA_F)*exp(ptilF)^-EPSILON_F+THETA_F*(exp(pF(-1))/exp(pF)*exp(pi(-1))^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi))^-EPSILON_F*exp(DelF(-1)); //E31
exp(mcF)=exp(rer)/exp(pF); //E32
exp(fF)=exp(ptilF)^-EPSILON_F*exp(yF)*exp(mcF)+BETTA/exp(a)^(SIGMA-1)*THETA_F*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi(+1)))^-EPSILON_F*(exp(ptilF)/exp(ptilF(+1)))^-EPSILON_F*(exp(pF)/exp(pF(+1)))^(-1-EPSILON_F)*exp(fF(+1))); //E33
exp(fF)=exp(ptilF)^(1-EPSILON_F)*exp(yF)*(EPSILON_F-1)/EPSILON_F+BETTA/exp(a)^(SIGMA-1)*THETA_F*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(pi)^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pi(+1)))^(1-EPSILON_F)*(exp(ptilF)/exp(ptilF(+1)))^(1-EPSILON_F)*(exp(pF)/exp(pF(+1)))^-EPSILON_F*exp(fF(+1))); //E34
exp(fW)=exp(mcW)*exp(wtil)^-EPSILON_W*exp(hd)+BETTA/exp(a)^(SIGMA-1)*THETA_W*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(gammaW)/exp(a)/exp(pi(+1)))^-EPSILON_W*(exp(wtil)/exp(wtil(+1)))^-EPSILON_W*(exp(w)/exp(w(+1)))^(-1-EPSILON_W)*exp(fW(+1))); //E35
exp(fW)=exp(wtil)^(1-EPSILON_W)*exp(hd)*(EPSILON_W-1)/EPSILON_W+BETTA/exp(a)^(SIGMA-1)*THETA_W*exp(varrho(+1))/exp(varrho)*(exp(lam(+1))/exp(lam)*(exp(gammaW)/exp(a)/exp(pi(+1)))^(1-EPSILON_W)*(exp(wtil)/exp(wtil(+1)))^(1-EPSILON_W)*(exp(w)/exp(w(+1)))^-EPSILON_W*exp(fW(+1))); //E36
1=(1-THETA_W)*exp(wtil)^(1-EPSILON_W)+THETA_W*(exp(w(-1))/exp(w)/exp(a(-1))*exp(gammaW(-1))/exp(pi))^(1-EPSILON_W); //E37
exp(DelW)=(1-THETA_W)*exp(wtil)^-EPSILON_W+THETA_W*(exp(w(-1))/exp(w)/exp(a(-1))*exp(gammaW(-1))/exp(pi))^-EPSILON_W*exp(DelW(-1)); //E38
exp(h)=exp(hd)*exp(DelW); //E39
exp(gammaW)=a_ss^ALPHA_W*exp(pi)^VARTHETA_W*pi_ss^(1-VARTHETA_W); //E40 or E48
exp(Theta)=exp(chitil)*(exp(c)-VARSIGMA*exp(c(-1))/exp(a(-1)))^-SIGMA; //E41 or E49
exp(chitil)=exp(chitil(-1))^(1-UPSILON)*(exp(c)-VARSIGMA*exp(c(-1))/exp(a(-1)))^(SIGMA*UPSILON); //E42 or E50

// Definitions 
exp(y) = exp(yNM) + exp(yCo);
exp(pYNM)*exp(yNM) = exp(pY)*exp(y)- exp(pF)*exp(pCostar)*exp(yCo);
exp(pCoyCo) = exp(rer)*exp(pCostar)*exp(yCo);
stb=tb/(exp(pY)*exp(y));
sdstar=exp(rer)*dstar/(exp(pY)*exp(y)); 
exp(gam_Y)=exp(pY)*exp(y)/(exp(pY(-1))*exp(y(-1)))*exp(a(-1));
exp(gam_YNM)=exp(pYNM)*exp(yNM)/(exp(pYNM(-1))*exp(yNM(-1)))*exp(a(-1));
exp(gam_C)=exp(c)/exp(c(-1))*exp(a(-1)); 
exp(gam_I)=exp(i)/exp(i(-1))*exp(a(-1)); 
exp(gam_W)=exp(w)/exp(w(-1))*exp(a(-1));  
exp(gam_K)=exp(k)/exp(k(-1))*exp(a(-1)); 
exp(gam_YC)=exp(yC)/exp(yC(-1))*exp(a(-1));
exp(gam_YF)=exp(yF)/exp(yF(-1))*exp(a(-1)); 
exp(gam_YH)=exp(yH)/exp(yH(-1))*exp(a(-1)); 
exp(gam_XF)=exp(xF)/exp(xF(-1))*exp(a(-1)); 
exp(gam_XH)=exp(xH)/exp(xH(-1))*exp(a(-1));
exp(gam_XHstar)=exp(xHstar)/exp(xHstar(-1))*exp(a(-1)); 
exp(gam_M)=exp(imp)/exp(imp(-1))*exp(a(-1)); 
exp(gam_LAM)=exp(lam)/exp(lam(-1))/exp(a(-1)); 
exp(piH)=exp(pH)/exp(pH(-1))*exp(pi); 
exp(piF)=exp(pF)/exp(pF(-1))*exp(pi); 

// Exogenous AR(1) processes (16)
varrho-log(varrho_ss)=RHO_varrho*(varrho(-1)-log(varrho_ss))+eps_varrho;
kappa-log(kappa_ss)=RHO_kappa*(kappa(-1)-log(kappa_ss))+eps_kappa; 
varomega-log(varomega_ss)=RHO_varomega*(varomega(-1)-log(varomega_ss))+eps_varomega;
z-log(z_ss)=RHO_z*(z(-1)-log(z_ss))+eps_z; 
a-log(a_ss)=RHO_a*(a(-1)-log(a_ss))+eps_a; 
yCo-log(yCo_ss)=RHO_yCo*(yCo(-1)-log(yCo_ss))+eps_yCo; 
g-log(g_ss)=RHO_g*(g(-1)-log(g_ss))+eps_g; 
zetao-log(zetao_ss)=RHO_zetao*(zetao(-1)-log(zetao_ss))+eps_zetao; 
zetau-log(zetau_ss)=RHO_zetau*(zetau(-1)-log(zetau_ss))+eps_zetau; 
Rstar-log(Rstar_ss)=RHO_Rstar*(Rstar(-1)-log(Rstar_ss))+eps_Rstar; 
ystar-log(ystar_ss)=RHO_ystar*(ystar(-1)-log(ystar_ss))+eps_ystar; 
pistar-log(pistar_ss)=RHO_pistar*(pistar(-1)-log(pistar_ss))+eps_pistar; 
pCostar-log(pCostar_ss)=RHO_pCostar*(pCostar(-1)-log(pCostar_ss))+eps_pCostar; 
eR-log(eR_ss)=RHO_eR*(eR(-1)-log(eR_ss))+eps_eR; 

// Observed variables
R_obs=100*(R-log(R_ss)); 
xi_obs=100*(xi-log(xi_ss)-(zetau-log(zetau_ss))); 
pi_obs=100*(pi-log(pi_ss)); 
rer_obs=100*(rer-log(rer_ss)); 
stb_obs=100*(stb-stb_ss); 
gam_Y_obs=100*(gam_Y-log(a_ss));
gam_YNM_obs=100*(gam_YNM-log(a_ss));
gam_C_obs=100*(gam_C-log(a_ss)); 
gam_I_obs=100*(gam_I-log(a_ss)); 
gam_W_obs=100*(gam_W-log(a_ss)); 
pCoyCo_obs=100*(pCoyCo-log(pCoyCo_ss));
yCo_obs=100*(yCo-log(yCo_ss)); 
g_obs=100*(g-log(g_ss)); 
Rstar_obs=100*(Rstar-log(Rstar_ss)); 
ystar_obs=100*(ystar-log(ystar_ss)); 
pistar_obs=100*(pistar-log(pistar_ss)); 
pCostar_obs=100*(pCostar-log(pCostar_ss)); 
hn_obs=100*(h-log(h_ss)); 
pi4_obs=100*(pi+pi(-1)+pi(-2)+pi(-3)-4*log(pi_ss));
piw4_obs=100*(gam_W+gam_W(-1)+gam_W(-2)+gam_W(-3)-4*log(a_ss))+pi4_obs;
piw_obs=100*(gam_W-log(a_ss))+pi_obs;
end;

%--------------------------------------------------------------------------
% 3. Steady state
%--------------------------------------------------------------------------
steady_state_model; 
// Computing the steady state and endogenous parameters
R_ss=a_ss^SIGMA*pi_ss/BETTA;
q_ss=1/varomega_ss;
rK_ss=q_ss*(a_ss^SIGMA/BETTA-1+DELTA);
piS_ss=a_ss^SIGMA*pi_ss/(BETTA*Rstar_ss*xi_ss);
pistar_ss=pi_ss/piS_ss;
ptilH_ss=1;
ptilF_ss=1;
gammaW_ss=a_ss^ALPHA_W*pi_ss;
wtil_ss=((1-THETA_W*(gammaW_ss/(a_ss*pi_ss))^(1-EPSILON_W))/(1-THETA_W))^(1/(1-EPSILON_W));
DelH_ss=ptilH_ss^-EPSILON_H;
DelF_ss=ptilF_ss^-EPSILON_F;
DelW_ss=((1-THETA_W)/(1-THETA_W*(gammaW_ss/(a_ss*pi_ss))^-EPSILON_W))*wtil_ss^-EPSILON_W;
hd_ss=h_ss/DelW_ss;
mcH_ss=(EPSILON_H-1)/EPSILON_H*ptilH_ss;
pm_ss=pH_ss*mcH_ss;
mcF_ss=(EPSILON_F-1)/EPSILON_F*ptilF_ss;
mcW_ss=((EPSILON_W-1)/EPSILON_W*(1-BETTA*a_ss^(1-SIGMA)*(gammaW_ss/(a_ss*pi_ss))^-EPSILON_W*THETA_W)/(1-BETTA*a_ss^(1-SIGMA)*(gammaW_ss/(a_ss*pi_ss))^(1-EPSILON_W)*THETA_W))*wtil_ss;
fW_ss=wtil_ss^-EPSILON_W*hd_ss*mcW_ss/(1-BETTA*a_ss^(1-SIGMA)*(gammaW_ss/(a_ss*pi_ss))^-EPSILON_W*THETA_W);
w_ss=(ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA)*pH_ss*mcH_ss*z_ss*a_ss^(1-ALPHA)/rK_ss^ALPHA)^(1/(1-ALPHA));
k_ss=a_ss*ALPHA*w_ss*hd_ss/((1-ALPHA)*rK_ss);
yH_ss=z_ss*(k_ss/a_ss)^ALPHA*(a_ss*hd_ss)^(1-ALPHA)/DelH_ss;
fH_ss=ptilH_ss^-EPSILON_H*yH_ss*mcH_ss/(1-BETTA*a_ss^(1-SIGMA)*THETA_H);
i_ss=k_ss*(1-(1-DELTA)/a_ss)/varomega_ss;
pF_ss=(1/O-(1-O)/O*pH_ss^(1-ETA))^(1/(1-ETA));
rer_ss=mcF_ss*pF_ss;
pYy_ss=pH_ss*yH_ss/(1-sCo_ss-pF_ss*(1-mcF_ss*DelF_ss)*O*pF_ss^-ETA*(1-stb_ss));
tb_ss=stb_ss*pYy_ss;
g_ss=sg_ss*pYy_ss;
yCo_ss=sCo_ss*pYy_ss/(rer_ss*pCostar_ss);
xHstar_ss=yH_ss-(1-O)*pH_ss^-ETA*(pYy_ss-tb_ss);
xH_ss=yH_ss-xHstar_ss;
xF_ss=(pH_ss*xHstar_ss+rer_ss*pCostar_ss*yCo_ss-tb_ss)/rer_ss;
yF_ss=xF_ss;
fF_ss=ptilF_ss^-EPSILON_F*yF_ss*mcF_ss/(1-BETTA*a_ss^(1-SIGMA)*THETA_F);
imp_ss=yF_ss*DelF_ss;
yC_ss=(xF_ss/O)*pF_ss^ETA;
c_ss=yC_ss-g_ss-i_ss;
y_ss=c_ss+i_ss+g_ss+xHstar_ss+yCo_ss-imp_ss;
pY_ss=(c_ss+i_ss+g_ss+tb_ss)/y_ss;
lam_ss=(c_ss-VARSIGMA*c_ss/a_ss)^(-SIGMA);
Theta_ss=1;
kappa_ss=mcW_ss*lam_ss*w_ss/(Theta_ss*h_ss^PHI); 
chitil_ss=c_ss^SIGMA*(1-VARSIGMA/a_ss)^SIGMA;
dstar_ss=-(tb_ss-(1-CHI)*rer_ss*pCostar_ss*yCo_ss)/(rer_ss*(1-Rstar_ss*xi_ss/pistar_ss/a_ss));
OSTAR=((1-VARKAPPA_X/a_ss)/(1-VARKAPPA_X))*(xHstar_ss/ystar_ss)*(pH_ss/rer_ss)^ETASTAR;
sdstar_ss=rer_ss*dstar_ss/pYy_ss;
pCoyCo_ss=  rer_ss*pCostar_ss*yCo_ss;
yNM_ss = y_ss -yCo_ss;
pYNM_ss = (pY_ss*y_ss - pF_ss*pCostar_ss*yCo_ss)/yNM_ss;

//Initial values for numerical solver
lam=log(lam_ss);
c=log(c_ss);
h=log(h_ss);
w=log(w_ss);
gammaW=log(gammaW_ss);
Theta=log(Theta_ss);
chitil=log(chitil_ss);
hd=log(hd_ss);
wtil=log(wtil_ss);
mcW=log(mcW_ss);
fW=log(fW_ss);
DelW=log(DelW_ss);
i=log(i_ss);
k=log(k_ss);
rK=log(rK_ss);
q=log(q_ss);
y=log(y_ss);
yNM=log(yNM_ss);
yC=log(yC_ss);
yF=log(yF_ss);
yH=log(yH_ss);
xF=log(xF_ss);
xH=log(xH_ss);
xHstar=log(xHstar_ss);
R=log(R_ss);
xi=log(xi_ss);
pi=log(pi_ss);
rer=log(rer_ss);
pH=log(pH_ss);
ptilH=log(ptilH_ss);
pF=log(pF_ss);
ptilF=log(ptilF_ss);
pY=log(pY_ss);
pYNM=log(pYNM_ss);
piS=log(piS_ss);
mcH=log(mcH_ss);
pm=log(pm_ss);
fH=log(fH_ss);
DelH=log(DelH_ss);
mcF=log(mcF_ss);
fF=log(fF_ss);
DelF=log(DelF_ss);
dstar=   (dstar_ss);
imp=log(imp_ss);
tb=   (tb_ss);
varrho=log(varrho_ss);
kappa=log(kappa_ss);
varomega=log(varomega_ss);
z=log(z_ss);
a=log(a_ss);
zetao=log(zetao_ss);
zetau=log(zetau_ss);
eR=log(eR_ss);
Rstar=log(Rstar_ss);
pistar=log(pistar_ss);
pCostar=log(pCostar_ss);
yCo=log(yCo_ss);
ystar=log(ystar_ss);
g=log(g_ss);
stb=   (stb_ss);
sdstar=   (sdstar_ss);
pCoyCo = log(pCoyCo_ss); 
gam_Y=log(a_ss);
gam_YNM=log(a_ss);
gam_C=log(a_ss);
gam_I=log(a_ss);
gam_W=log(a_ss);  
gam_K=log(a_ss);
gam_YC=log(a_ss);
gam_YF=log(a_ss);
gam_YH=log(a_ss);
gam_XF=log(a_ss);
gam_XH=log(a_ss);
gam_XHstar=log(a_ss);
gam_M=log(a_ss);
gam_LAM=log(1/a_ss);
piH=log(pi_ss);
piF=log(pi_ss);
gam_Y_obs=0;
gam_YNM_obs=0;
gam_C_obs=0;
gam_I_obs=0;
gam_W_obs=0;
R_obs=0;
xi_obs=0;
pi_obs=0;
rer_obs=0;
yCo_obs=0;
pCoyCo_obs=0;
g_obs=0;
Rstar_obs=0;
ystar_obs=0;
pistar_obs=0;
pCostar_obs=0;
stb_obs=0;
hn_obs=0;
end;

//Call steady state solver
steady;

//Check Blanchard-Kahn conditions
check;

%--------------------------------------------------------------------------
% 4. Estimation and simulation
%--------------------------------------------------------------------------
// Observed variables
varobs gam_YNM_obs gam_C_obs gam_I_obs gam_W_obs R_obs xi_obs pi_obs rer_obs
pCoyCo_obs g_obs Rstar_obs pistar_obs ystar_obs pCostar_obs hn_obs;

shocks;
// Structural shocks
var eps_varrho; stderr SIG_varrho;
var eps_kappa;  stderr SIG_kappa;
var eps_varomega; stderr SIG_varomega;
var eps_z;        stderr SIG_z;
var eps_a; 		  stderr SIG_a;
var eps_zetao;    stderr SIG_zetao;
var eps_zetau;    stderr SIG_zetau;
var eps_eR;       stderr SIG_eR;
var eps_yCo;      stderr SIG_yCo;
var eps_g;        stderr SIG_g;
var eps_Rstar;    stderr SIG_Rstar;
var eps_ystar;    stderr SIG_ystar;
var eps_pistar;   stderr SIG_pistar;
var eps_pCostar;  stderr SIG_pCostar;
var eps_aux;      stderr eps;
var eps_aux1;     stderr eps;

// Measurement errors, calibrated 
var pCoyCo_obs; stderr sqrt(41.3477^2*0.1);   
var gam_YNM_obs; 	stderr sqrt(1.1437^2*0.1);   
var gam_C_obs; 	stderr sqrt(1.2061^2*0.1);
var gam_I_obs; 	stderr sqrt(3.7542^2*0.1);
var gam_W_obs; 	stderr sqrt(0.5847^2*0.1);      
var R_obs; 		stderr sqrt(0.4164^2*0.0);
var xi_obs; 	stderr sqrt(0.1467^2*0.1);
var pi_obs; 	stderr sqrt(0.6916^2*0.1);
var rer_obs; 	stderr sqrt(5.1679^2*0.1); 
var hn_obs; 	stderr sqrt(1.8690^2*0.1);
end;

estimated_params;
// Deep parameters     
	UPSILON,    beta_pdf, 0.5, 0.25;
    PHI, , 0, , normal_pdf, 2, 2; 
    VARSIGMA,   beta_pdf, 0.7, 0.1;
    PSI,        inv_gamma_pdf, 0.01, inf;
    ETA,        inv_gamma_pdf, 1.5, 0.25;
    ETASTAR,    inv_gamma_pdf, 0.25, 0.1;
    GAMA,       normal_pdf, 4, 1.5; 
	THETA_W,    beta_pdf, 0.5, 0.1;
	VARTHETA_W, beta_pdf, 0.5, 0.15;
    THETA_H,    beta_pdf, 0.5, 0.1;
    VARTHETA_H, beta_pdf, 0.5, 0.15;
    THETA_F,    beta_pdf, 0.5, 0.1;
    VARTHETA_F, beta_pdf, 0.5, 0.15;
    RHO_R,      beta_pdf, 0.75, 0.1;
    ALPHA_PI, , 1, , normal_pdf, 1.5, 0.1;
    ALPHA_Y,    normal_pdf, 0.5/4, 0.05;
// AR coefficients
    RHO_yCo,    beta_pdf, 0.75, 0.1;
    RHO_varrho, beta_pdf, 0.75, 0.1;
    RHO_kappa,  beta_pdf, 0.75, 0.1;
    RHO_varomega, beta_pdf, 0.75, 0.1;
    RHO_z,        beta_pdf, 0.75, 0.1;
    RHO_a,        beta_pdf, 0.75/2, 0.1;
    RHO_zetao,    beta_pdf, 0.75, 0.1;
    RHO_zetau,    beta_pdf, 0.75, 0.1;
// Innovation std. deviations
    stderr eps_yCo,    inv_gamma_pdf, 0.01, inf;
    stderr eps_varrho, inv_gamma_pdf, 0.01, inf;
    stderr eps_kappa,  inv_gamma_pdf, 0.01, inf;
    stderr eps_varomega, inv_gamma_pdf, 0.01, inf;
    stderr eps_z,     inv_gamma_pdf, 0.01, inf;
    stderr eps_a,     inv_gamma_pdf, 0.01, inf;
    stderr eps_zetao, inv_gamma_pdf, 0.01/4, inf;
    stderr eps_zetau, inv_gamma_pdf, 0.01/4, inf;
    stderr eps_eR,    inv_gamma_pdf, 0.01/4, inf;
end;


//The following computes the posterior mode
//estimation(datafile=estim_data, xls_sheet=estim_data, mode_compute=4, mh_replic= 0 plot_priors=1); 

//The following runs MCMC replications to estimate posterior distribution. It uses the posterior mode as initial point
//estimation(datafile=estim_data, xls_sheet=estim_data, mode_compute=0, mode_file=standar_mode, mh_replic=0, mh_jscale=0.35, mh_nblocks=1, plot_priors=0);

//The following simulates the model by first-order perturbation
stoch_simul(order=1,periods=0,irf=40,nograph);

//The following computes recursive forecasts
/*
@#for n in 22:59
    estimation(datafile=estim_data,xls_sheet=estim_data,plot_priors=0,first_obs=1,nobs=@{n}, mh_replic=0,forecast=10,nograph) gam_YNM_obs pi_obs R_obs rer_obs hn_obs gam_W_obs;
    save Standar_@{n}.mat       
@#endfor
*/