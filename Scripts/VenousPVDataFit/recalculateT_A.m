clear; figure(1); clf

 
% Parameters

alpha =11.5; % unitless

beta = 1.5; % unitless

a = 159.9; % Pa

b = 7.99e-3; % Pa

c = 3.49e3; % Pa

Vn = 1e-4; % m^3 100 mL (set to arbitrary value)

gamma = 0.40; % unitless (set to arbitrary value)

len = 0.1; % m (length of vessel)

V0 = gamma*Vn;

L0 = sqrt(4*pi*V0/len); % calculate L0 from V0

 

V = (1:0.1:4)*V0; % range of volumes (m^3)

 

% Computing pressure at each V(i), with phi 0

for i = 1:length(V)

  L = sqrt(4*pi*V(i)/len); % circumference

  r = L/(2*pi); % radius

  f = L*(L-L0)/L0;

  g = L0*(exp(alpha*((L-L0))/L0)-1);

  h = (L-L0)*exp(-beta*((L-L0)/L0 )^2);

  T_P = a*f + b*g;

  T_A = max(0,(0-0.2)/0.8)*c*h;

  P(i) = (T_P + T_A)/r;

end

 

% plot volume (in mL) vs. pressure (mmHg)

figure(1); plot(V.*1e6,P*0.00750062); grid on; hold on;

ylabel('Pressure (mmHg)');

xlabel('Volume (mL)');

 

% Computing pressure at each V(i), with phi 0.25

for i = 1:length(V)

  L = sqrt(4*pi*V(i)/len); % circumference

  r = L/(2*pi); % radius

  f = L*(L-L0)/L0;

  g = L0*(exp(alpha*((L-L0))/L0)-1);

  h = (L-L0)*exp(-beta*((L-L0)/L0 )^2);

  T_P = a*f + b*g;

  T_A = max(0,(0.25-0.2)/0.8)*c*h;

  P(i) = (T_P + T_A)/r;

end

 

% plot volume (in mL) vs. pressure (mmHg)

figure(1); plot(V.*1e6,P*0.00750062); grid on; hold on;

ylabel('Pressure (mmHg)');

xlabel('Volume (mL)');

 

% Computing pressure at each V(i), with phi 1.0

for i = 1:length(V)

  L = sqrt(4*pi*V(i)/len); % circumference

  r = L/(2*pi); % radius

  f = L*(L-L0)/L0;

  g = L0*(exp(alpha*((L-L0))/L0)-1);

  h = (L-L0)*exp(-beta*((L-L0)/L0 )^2);

  T_P = a*f + b*g;

  T_A = max(0,(1.0-0.2)/0.8)*c*h;

  P(i) = (T_P + T_A)/r;

end

 

% plot volume (in mL) vs. pressure (mmHg)

figure(1); plot(V.*1e6,P*0.00750062); grid on; hold on;

ylabel('Pressure (mmHg)');

xlabel('Volume (mL)');