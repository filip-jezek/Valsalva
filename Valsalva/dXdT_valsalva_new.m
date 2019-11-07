function [dxdt] = dXdT_valsalva_new(t,x,pars)


%% Parameters 

ts    = pars(1); 
te    = pars(2); 
K_Pth  = pars(3); 
q_Pth = pars(4); 

phi1  = pars(5); 
K_phi1 = pars(6); 
K_phi2 = pars(7); 
q_phi1 = pars(8); 
q_phi2 = pars(9); 

alphaR = pars(10);
alphaC = pars(11);

R_ta  = pars(12);
R_sa1 = pars(13);
R_sv  = pars(14);
R_tv  = pars(15);
R_pa  = pars(16);
R_pv  = pars(17);
Rv    = pars(18);

C_ta  = pars(19);
C_sa  = pars(20);
C_sv1 = pars(21);
C_tv  = pars(22);
C_pa  = pars(23);
C_pv  = pars(24);

ep = 1;

%% State variables

theta = x(1); % beat counter
V_ta  = x(2); % volume thoracic aorta (mL)
V_sa  = x(3); % volume systemic arteries, outside of TC (mL)
V_sv  = x(4); % volume systemic veins, outside of TC (mL)
V_tv  = x(5); % volume thoracic vena cava (mL)
V_rv  = x(6); % volume right ventricle (mL)
V_pa  = x(7); % volume pulmonary arteries (mL)
V_pv  = x(8); % volume pulmonary veins (mL)
V_lv  = x(9); % volume thoracic LV (mL)

%% Auxiliary Equations

% Thoracic pressure and phi
if t < ts
  Pth = 0;
  phi = phi1; % sympathetic tone
else
  if t < te
    Pth = K_Pth*(1 - exp( -q_Pth*(t-ts) ));
    phi = phi1 + K_phi1*(1 - exp( -q_phi1*(t-ts) )); % sympathetic tone
  else
    Pth = K_Pth*exp( -q_Pth*(t-te) );
    phi = phi1 + K_phi2*exp( -q_phi2*(t-te) );
  end
end

H = 1 + 1*(phi-phi1);

R_sa = R_sa1 * (1 + alphaR*(phi-phi1));
C_sv = C_sv1 / (1 + alphaC*(phi-phi1));

% Pressures
P_ta = V_ta/C_ta + Pth; % thoracic aorta
P_sa = V_sa/C_sa; % systemic arteries
P_sv = V_sv/C_sv; % systemic veins
P_tv = V_tv/C_tv + Pth; % thoracic veins
P_pa = V_pa/C_pa + Pth; % pul. artery
P_pv = V_pv/C_pv + 1.05*Pth; % pul. vein
P_lv = PVfunction_new(theta,V_lv,ep,pars) + Pth;
P_rv = PVfunction_new(theta,V_rv,0.3,pars) + Pth;

F_lv = max(0, (P_lv-P_ta)/Rv);
F_ta = (P_ta-P_sa)/R_ta;
F_sa = (P_sa-P_sv)/R_sa; 
F_sv = (P_sv-P_tv)/R_sv; 
F_tv  = max(0, (P_tv-P_rv)/R_tv);
F_rv = max(0, (P_rv-P_pa)/Rv);
F_pa = (P_pa-P_pv)/R_pa; 
F_pv  = max(0, (P_pv-P_lv)/R_pv);

%% DEs

dtheta = H;
dV_ta = F_lv - F_ta; % V_ta
dV_sa = F_ta - F_sa; % V_sa
dV_sv = F_sa - F_sv; % V_sv
dV_tv = F_sv - F_tv; %V_tv
dV_rv = F_tv - F_rv; %V_rv
dV_pa = F_rv - F_pa; %V_pa
dV_pv = F_pa - F_pv; %V_pv
dV_lv = F_pv - F_lv; %V_lv

dxdt = [dtheta; dV_ta; dV_sa; dV_sv; 
    dV_tv; dV_rv; dV_pa; dV_pv; dV_lv]; 
