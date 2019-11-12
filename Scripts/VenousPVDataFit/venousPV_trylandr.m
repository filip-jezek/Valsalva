
%Given from Cell ML model
l  = 45; %cm 
rn = .001; %cm 

%Assumptions 
gamma = 0.5; 
phibar = 0.25; 

%Nominal values 
Ln = 2*pi*rn; 
Vn = pi*rn^2*l; 

%Zero pressure values
V0 = gamma*Vn;
L0 = 2*pi*sqrt(V0/pi/l); 

%Max volume 
VM = 4*V0; 
LM = 2*pi*sqrt(VM/pi/l);

%Functional components 
f = @(L) L.*(L - L0)/L0; 
g = @(L) L0*(exp( 11.5*(L - L0)/L0 ) - 1); 
h = @(L) (L - L0).*exp( -1.5*((L - L0)/L0).^2 ); 

a = 1.2; 
b = 6e-5; 
c = 26.25; 

T_P      = @(L) a*f(L) + b*g(L);
T_A      = @(L) c*h(L); 
T_nom    = @(L) T_P(L) + phibar*T_A(L);
T_maxact = @(L) T_P(L) + T_A(L);

P_nom    = @(L) T_nom(L)*(2*pi)./L;
P_pas    = @(L) T_P(L)*(2*pi)./L; 
P_act    = @(L) T_A(L)*(2*pi)./L; 
P_maxact = @(L) T_maxact(L).*(2*pi)./L; 

P_nom(Ln) 
P_pas(Ln)
P_act(Ln)
P_maxact(Ln)

P_nom(LM) 
P_pas(LM)
P_act(LM)
P_maxact(LM)

P_nom(L0) 
P_pas(L0)
P_act(L0)
P_maxact(L0)
