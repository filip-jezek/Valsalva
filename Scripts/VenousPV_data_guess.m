% assumptions
% Vn - known at nominal activation
Vn = 1.25*V0; % at passive
Vm = 3*V0 % at passive
%{
Know - Pn,Vn at phi = 0.25 
Assume - 
1. At phi = 0: P0 = 0,V0 = gamma*Vn for gamma = 0.5
2. At phi = 1: PM = 30/40, VM = 3.5*V0 = 3.5*gamma*Vn 
3. T_phi0 *alpha = T_phi1 for alpha = 2.5 at Vn
Maybe 4. Choose offset value for Gaussian 
max activation tension must be 0 at diameter 0
T_A = c_A * (L/Lstar) * exp( -k ((L - Lstar)/Lstar)^2)
Vstar = beta*Vn,for Vstar the peak max active and assume beta = 1 for the time being

Find k using expression for TA only and assign Lstar and get a value for c_A but don't use it 
 Come up with algebraic expressions for scaling values for TP and TA
T_P + phi*T_A @ => fit k, c_A

%}

%% fit of T_A from Vanhoute
recalculate_PV_data_into_tensions;
%% 
length = 10e-3; % the paper is not specific about lengths of the venous segments. Lets have a guess here [m]
gamma = 0.5; % V0 = gamma*Vn
beta = 0.5; % Vn = beta*Vstar

% resulting fit: x_data_L to act_tension_data
% General model:
%      f(L) = c_A * (L/Lstar) * exp( -k* ((L - Lstar)/Lstar)^2)
% Coefficients (with 95% confidence bounds):
%        Lstar =     0.01819  (0.01725, 0.01913)
%        c_A =       15.03  (13.68, 16.38)
%        k =       3.425  (2.148, 4.701)
% 
% Goodness of fit:
%   SSE: 10.58
%   R-square: 0.9479
%   Adjusted R-square: 0.9363
%   RMSE: 1.084

%  T_A = @(L) c_A * (L/Lstar) .* exp( -k* ((L - Lstar)/Lstar).^2);
% get rid of the linear term - although the tension is not going to be zero
% at zero circumferential lengths, it shouldnt cause much troubles.
% If it did - multiply with something along (1-(exp(-L/Lstar))) 
% General model:
%      f(L) = c_A * exp( -k* ((L - Lstar)/Lstar)^2)
% Coefficients (with 95% confidence bounds):
%        Lstar =     0.02062  (0.01999, 0.02124)
%        c_A =       16.05  (14.89, 17.21)
%        k =       5.002  (3.771, 6.233)

% fit to act_tension_data;
       Lstar =     0.02062 ;
       c_A =       16.05  ;
       k =       5.002  ;

T_A = @(L) c_A * exp( -k* ((L - Lstar)/Lstar).^2); 
% backcalculate Vn for estimation of Lstar
rstar = Lstar/2/pi;
Vstar = pi*length*rstar^2;
Vn = Vstar*beta;
V0 = Vn*gamma;
Ln = 2*pi*sqrt(Vn/length/pi);
L0 = 2*pi*sqrt(V0/length/pi);

y_bnds = [ min(act_tension_data), max(act_tension_data)*1.1];
figure(201); clf;hold on;
plot(x_data_L, [0;act_tension_data(2:end)], 'r-*');
plot(x_samples_L, T_A(x_samples_L), 'b-');
plot([1 1]*Lstar, y_bnds,'b--', 'LineWidth', 2);
plot([1 1]*Ln, y_bnds, 'm--');
plot([1 1]*L0, y_bnds,'k:');
legend('data', 'fit', 'L0', 'Ln', 'Lstar');
xlabel('Circumferential length [m]')
ylabel('Tension [N]')

figure(202); clf;hold on;
y_bnds = [0, 40];
plot(x_data_L, [0;act_tension_data(2:end)]./x_data_r/133.32, 'r-*');
plot(x_samples_L, T_A(x_samples_L)./x_samples_r/133.32 , 'b-');
plot(x_samples_L, T_A(x_samples_L)./x_samples_r/133.32*0.2 , 'b--');
plot([1 1]*Lstar, y_bnds,'b--', 'LineWidth', 2);
plot([1 1]*Ln, y_bnds, 'm--');
plot([1 1]*L0, y_bnds,'k:');
legend('data', 'fit', 'baseline', 'Lstar', 'Ln', 'Lstar');
xlabel('Circumferential length [m]')
ylabel('Pressure [mmHg]')