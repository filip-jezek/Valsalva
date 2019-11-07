%% Pressure/Volume function ventricle

function [P_lv] = PVfunction_new(Theta,V_lv,ep,pars)

  Vd   = 0;    % unstressed volume (ml)
  
  E_max = pars(25); 
  E_min = pars(26); 
  T_Mf  = pars(27);
  T_Rf  = pars(28); 
  
  EM = ep*E_max;
  
  tTilde = mod(Theta,1); % fraction of cardiac cycle

  if tTilde < T_Mf
    E = (EM - E_min)*(1 - cos(pi*tTilde/T_Mf))/2 + E_min;
  elseif (tTilde >= T_Mf) && (tTilde < (T_Mf + T_Rf))
    E = (EM - E_min)*(cos(pi*(tTilde - T_Mf)/T_Rf) + 1)/2 + E_min;
  else 
    E = E_min;
  end

P_lv = E*( V_lv - Vd );
% P_lv = E*Vd*( (V_lv/Vd-1) + 0.2*(V_lv/Vd-1)^2 );