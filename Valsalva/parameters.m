function [pars, x0] = parameters

CO = 5000/60; 

%% Initial conditions
theta = 0;
V_ta0 = 5;
V_sa0 = 75;
V_sv0 = 3000;
V_tv0 = 465;
V_rv0 = 120;
V_pa0 = 100;
V_pv0 = 671;
V_lv0 = 125;

x0 = [theta V_ta0 V_sa0 V_sv0 V_tv0 V_rv0 V_pa0 V_pv0 V_lv0]';

%% Parameters
ts    = 30; 
te    = 49; 
K_Pth = 40; 
q_Pth = 2; 

phi1   = 0.25; 
K_phi1 = 0.5; 
K_phi2 = 0.4; 
q_phi1 = 0.1; 
q_phi2 = 0.5; %or 0.1? 

alphaR = 1.5;
alphaC = 6;

R_ta  = (100-98)/CO;
R_sa1 = (98-5)/CO;
R_sv  = (5-1)/CO;
R_tv  = 1/CO;
R_pa  = (18-1)/CO;
R_pv  = 1/CO;
Rv    = 0.001;

C_ta  = 7/100;
C_sa  = 100/98 ;
C_sv1 = 3000/10; 
C_tv  = 75/1;
C_pa  = 100/18;
C_pv  = 100/1;

E_max = 2.5; %2.6; 
E_min = 0.05; %0.008;
T_Mf  = 0.3;    % time to max E
T_Rf  = 0.15;   % relaxation time

pars = [ts, te, K_Pth, q_Pth, ... 
    phi1, K_phi1, K_phi2, q_phi1, q_phi2, ...
    alphaR, alphaC, ...
    R_ta, R_sa1, R_sv, R_tv, R_pa,R_pv, Rv, ...
    C_ta, C_sa, C_sv1, C_tv, C_pa, C_pv, ...
    E_max, E_min, T_Mf, T_Rf]'; 