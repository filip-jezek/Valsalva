clear

%% Parameters

[pars,x0] = parameters; 

ts    = pars(1); 
te    = pars(2); 
K_Pth = pars(3); 
q_Pth = pars(4); 

phi1   = pars(5);
K_phi1 = pars(6);
K_phi2 = pars(7);
q_phi1 = pars(8);
q_phi2 = pars(9); 

alphaR = pars(10);
alphaC = pars(11);

C_ta  = pars(19);
C_sa  = pars(20);
C_sv1 = pars(21);
C_tv  = pars(22);
C_pa  = pars(23);
C_pv  = pars(24);

%% Simulation

options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[t,x] = ode15s(@dXdT_valsalva_new,[0 80],x0,options,pars);

theta = x(:,1);
V_ta  = x(:,2); % volume thoracic aorta (mL)
V_sa  = x(:,3); % volume systemic arteries, outside of TC (mL)
V_sv  = x(:,4); % volume systemic veins, outside of TC (mL)
V_tv  = x(:,5); % volume thoracic vena cava (mL)
V_rv  = x(:,6); % volume right ventricle (mL)
V_pa  = x(:,7); % volume pulmonary arteries (mL)
V_pv  = x(:,8); % volume pulmonary veins (mL)
V_lv  = x(:,9); % volume thoracic LV (mL)
V_total = V_ta + V_sa + V_sv + V_tv + V_rv + V_pa + V_pv + V_lv;

%% Equations

Pth = zeros(size(t)) + ...
    ((t>ts).*(t<te)).*(K_Pth*(1 - exp( -q_Pth*(t-ts) ))) + ...
    (t>=te).*(K_Pth*(exp( -q_Pth*(t-te) )));

phi = phi1*ones(size(t)) + ...
    ((t>ts).*(t<te)).*(K_phi1*(1 - exp( -q_phi1*(t-ts) ))) + ...
    (t>=te).*(K_phi2*(exp( -q_phi2*(t-te) ))); 

C_sv = C_sv1 ./ (1 + alphaC*(phi-phi1));

% Pressures
P_ta = V_ta./C_ta + Pth;
P_sa = V_sa./C_sa;
P_sv = V_sv./C_sv;
P_tv = V_tv./C_tv + Pth;
P_pa = V_pa./C_pa + Pth;
P_pv = V_pv./C_pv + 1.05*Pth;

for i = 1:length(t)
  P_lv(i) = PVfunction_new(theta(i),V_lv(i),1,pars) + Pth(i);
  P_rv(i) = PVfunction_new(theta(i),V_rv(i),0.3,pars) + Pth(i);
end

%% Plots

figure(1)
plot(t,P_sa,t,P_pa)
legend('P_{sa}','P_{pa}');

print -dpng Figures/psa_new.png

figure(2)
plot(t,V_lv,t,V_rv)
legend('V_{lv}','V_{rv}');

print -dpng Figures/vlv_new.png

figure(3); plot(t,P_sv,t,P_pv)
legend('P_{sv}','P_{pv}');

print -dpng Figures/psv_new.png 

figure(4); plot(t,P_lv,t,P_rv)
legend('P_{lv}','P_{rv}');

print -dpng Figures/plv_new.png

figure(5)
plot(t,60 + 60*(phi-phi1))
axis([0 80 0 120])
legend('HR')

print -dpng Figures/hr_new.png