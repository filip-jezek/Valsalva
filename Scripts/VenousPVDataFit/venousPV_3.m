
clear all
plotAll = false;

%% Data
%{
Load and plot the data
%}

load Moreno70_data.mat 
load Vanhoutte69_data.mat 

if plotAll
    % Moreno
    figure(1) 
    clf
    hold on 
    plot([min(Pd_Moreno) max(Pd_Moreno)],zeros(2,1),'k')
    plot(zeros(2,1),[min(Vd_sca_Moreno) max(Vd_sca_Moreno)],'k')
    plot(Pd_Moreno,Vd_sca_Moreno,'mo','MarkerSize',5,'linewidth',2)
    xlabel('Pressure (mmHg)')
    ylabel('Relative volume')
    set(gca,'FontSize',15)
    title('Moreno data - original')

    % Vanhoutte
    figure(2)
    clf
    hold on 
    plot([min(Vd_pass_Vanh) max(Vd_pass_Vanh)],zeros(2,1),'k')
    plot(zeros(2,1),[min(Pd_pass_Vanh) max(Pd_maxact_Vanh)],'k')
    plot(Vd_pass_Vanh,Pd_pass_Vanh,'bo','MarkerSize',5,'linewidth',2)
    plot(Vd_maxact_Vanh,Pd_maxact_Vanh,'ko','color',[.5 0 .5],'MarkerSize',5,'linewidth',2)
    plot(Vd_pass_Vanh,Pd_maxact_Vanh - Pd_pass_Vanh,'ro','MarkerSize',5,'linewidth',2)
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',15)
    title('Vanhoutte data - original')
end

%% Parameters and target values

phibar = 0.25; % normal autonomic outflow 

% phi to activation conversion
veins_phi_no = 0.2;
phi_ns = 1/((1 - veins_phi_no));
A = min(max(0.001, (phibar - veins_phi_no)*phi_ns), 1)

P_n = 4.4; %(max(Pd_Moreno) - min(Pd_Moreno)) * A;    % mmHg
V_n = 50;   % mL - From Filip's simulations' approximations 
l   = 20;   % cm - From Cell ML model 

gamma = 0.5;  % scaling factor
V_0   = gamma * V_n; % mL

%% Transform Moreno data using V_0

Vd_Moreno = V_0 * Vd_sca_Moreno + V_0; % rescale data based on prescribed V_0 (m^3)

if plotAll
    figure(3) 
    clf
    hold on 
    plot([min(Pd_Moreno) max(Pd_Moreno)],zeros(2,1),'k')
    plot([min(Pd_Moreno) max(Pd_Moreno)],V_0 * ones(2,1),'k:')
    plot([min(Pd_Moreno) max(Pd_Moreno)],V_n * ones(2,1),'k:')
    plot(zeros(2,1),[min(Vd_Moreno) max(Vd_Moreno)],'k')
    plot(P_n * ones(2,1),[min(Vd_Moreno) max(Vd_Moreno)],'k:')
    plot(Pd_Moreno,Vd_Moreno,'mo','MarkerSize',5,'linewidth',2)
    xlabel('Pressure (mmHg)')
    ylabel('Volume (mL)')
    set(gca,'FontSize',15)
    title('Moreno data - volume scaled')
end

%% Transform Vanhoutte volume data to be on the same scale as the transformed Moreno data

Pd_act_Vanh = Pd_maxact_Vanh - Pd_pass_Vanh; 

Amp_Moreno = max(Vd_Moreno) - V_0; 
Amp_Vanh   = max(Vd_pass_Vanh) - min(Vd_pass_Vanh); 

Vd_pass_sca = Amp_Moreno/Amp_Vanh * Vd_pass_Vanh + V_0; 
Vd_maxact_sca = Amp_Moreno/Amp_Vanh * Vd_maxact_Vanh + V_0; 
Vd_act_sca = Vd_pass_sca; 

if plotAll,
    figure(4)
    clf
    hold on 
    plot([min(Pd_Moreno) max(Pd_Moreno)],zeros(2,1),'k')
    plot(zeros(2,1),[min(Vd_Moreno) max(Vd_Moreno)],'k')
    plot(Pd_Moreno,Vd_Moreno,'mo','MarkerSize',5,'linewidth',2)
    plot(Pd_pass_Vanh,Vd_pass_sca,'bo','MarkerSize',5,'linewidth',2)
    title('Scaled Moreno data and scaled passive Vanhoutte data')
    set(gca,'FontSize',15)
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')

    figure(5)
    clf
    hold on 
    plot([min(Vd_pass_sca) max(Vd_pass_sca)],zeros(2,1),'k')
    plot(zeros(2,1),[min(Pd_pass_Vanh) max(Pd_maxact_Vanh)],'k')
    plot(Vd_pass_sca,Pd_pass_Vanh,'bo','MarkerSize',5,'linewidth',2)
    plot(Vd_maxact_sca,Pd_maxact_Vanh,'ko','color',[.5 0 .5],'MarkerSize',5,'linewidth',2)
    plot(Vd_pass_sca,Pd_act_Vanh,'ro','MarkerSize',5,'linewidth',2)
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',15)
    title('Vanhoutte data - volume scaled') 
end


%% Convert scaled Moreno data to length-tension 

% convert pressure from mmHg to Pa
mmHg2SI = 133.322;
Pd_Moreno = Pd_Moreno * mmHg2SI;

% convert V_0 from mL to m^3
V_0 = V_0 * 1e-6;
Vd_Moreno = Vd_Moreno * 1e-6; 

% convert cm to m
l = l * 1e-2;  

L_0 = 2 * pi * sqrt(V_0 / (l * pi)); 

L_n = 2 * pi * sqrt( V_n * 1e-6 / (l * pi)); 
Ld_Moreno = 2 * pi * sqrt( Vd_Moreno / (l * pi)); % transform volume to circumferential length (m)  
Td_Moreno = Pd_Moreno .* Ld_Moreno / (2 * pi); % transform pressure and circumferential length to tension via Law of Laplace (

figure(6) 
clf
hold on 
plot([min(Ld_Moreno) max(Ld_Moreno)],zeros(2,1),'k')
h1 = plot(Ld_Moreno, Td_Moreno,'ko','MarkerSize',5,'linewidth',2);
xlabel('Length (m)')
ylabel('Tension (Pa)')
set(gca,'FontSize',15)

%{
Used cftool on Ld_Moreno and Td_Moreno to fit parameter d = 15
(dimensionless) 
%}
d = 11.5; 

%% Convert scaled Vanhoutte data to length-tension 

% convert pressure from mmHg to Pa 
Pd_pass_Vanh = Pd_pass_Vanh * mmHg2SI ; 
Pd_maxact_Vanh = Pd_maxact_Vanh * mmHg2SI; 
Pd_act_Vanh = Pd_act_Vanh * mmHg2SI; 

% convert volumes from mL to m^3 
Vd_pass_sca = Vd_pass_sca * 1e-6; 
Vd_maxact_sca = Vd_maxact_sca * 1e-6; 
Vd_act_sca = Vd_act_sca * 1e-6; 

Ld_pass_Vanh   = 2 * pi * sqrt( Vd_pass_sca / (l * pi)); 
Ld_maxact_Vanh = 2 * pi * sqrt( Vd_maxact_sca / (l * pi)); 
Ld_act_Vanh    = 2 * pi * sqrt( Vd_act_sca / (l * pi)); 

Td_pass_Vanh   = Pd_pass_Vanh   .* Ld_pass_Vanh   / (2 * pi); 
Td_maxact_Vanh = Pd_maxact_Vanh .* Ld_maxact_Vanh / (2 * pi); 
Td_act_Vanh    = Pd_act_Vanh    .* Ld_act_Vanh    / (2 * pi); 

figure(6)
hold on 
plot(zeros(2,1),[min(Td_Moreno) max(Td_maxact_Vanh)],'k')
plot(L_0 * ones(2,1),[min(Td_Moreno) max(Td_maxact_Vanh)],'k:')
plot(L_n * ones(2,1),[min(Td_Moreno) max(Td_maxact_Vanh)],'k:')
h2 = plot(Ld_pass_Vanh, Td_pass_Vanh,'bo','MarkerSize',5,'linewidth',2);
h3 = plot(Ld_maxact_Vanh, Td_maxact_Vanh,'ko','color',[.5 0 .5],'MarkerSize',5,'linewidth',2);
h4 = plot(Ld_act_Vanh, Td_act_Vanh,'ro','MarkerSize',5,'linewidth',2);
xlabel('Length (m)')
ylabel('Tension (Pa)')
set(gca,'FontSize',15) 

%{
Used cftool on Ld_act_Vanh and Td_act_Vanh to fit parameter e = 1.5
(dimensionless)
%} 
e = 1;  

%% Model 

% Tension ratio T_maxact/T_pass at L_n 
[~,i_pass] = min(abs(L_n - Ld_pass_Vanh));
[~,i_maxact] = min(abs(L_n - Ld_maxact_Vanh));
alpha = Td_maxact_Vanh(i_maxact) / Td_pass_Vanh(i_pass);
alpha = round(alpha); 

L = 0:.001:Ld_pass_Vanh(end); 

T_n = P_n * mmHg2SI * L_n / (2 * pi);  

P_M = max(Pd_maxact_Vanh/ mmHg2SI ); 
V_M = 4 * V_0; 
L_M = 2 * pi * sqrt(V_M / (l * pi)); 
T_M = P_M * mmHg2SI * L_M / (2 * pi); 

f = @(L) L .* (L - L_0) / L_0; 
g = @(L) L .* (exp( d * (L - L_0) ./ L_0) - 1); 
h = @(L) (L - L_0) .* exp( -e * ((L - L_0)./L_0).^2); 

fhat = f(L_M)/f(L_n); 
ghat = g(L_M)/g(L_n); 
hhat = h(L_M)/h(L_n);   

a = (- T_M + T_n / (1 + A * (alpha - 1)) * (ghat + (alpha - 1) * hhat)) ...
    / (f(L_n) * (ghat - fhat));
b = (T_M - T_n / (1 + A * (alpha - 1)) * (fhat + (alpha - 1) * hhat)) ... 
    / (g(L_n) * (ghat - fhat));
c = T_n / (1 + A * (alpha - 1)) * (alpha - 1) / h(L_n); 

sprintf("a %d, b %d, c %d", a, b, c)

%%
T_P  = @(L) a * f(L) + b * g(L); 
T_A  = @(L) c * h(L); 

% build the characteristics
for i = 1:length(L)
    if L(i) < L_0 
        Tpas(i) = T_P(L(i));  
        Tact(i) = T_P(L(i));  
        Tbar(i) = T_P(L(i));  
        Tmax(i) = T_P(L(i));  
    else 
        Tpas(i) = T_P(L(i)) + 0 * T_A(L(i)); 
        Tact(i) = T_A(L(i));  
        Tbar(i) = T_P(L(i)) + A * T_A(L(i)); 
        Tmax(i) = T_P(L(i)) + 1 * T_A(L(i)); 
    end 
end 
    
figure(6)
hold on 
plot(L,Tpas,'b','linewidth',2)
plot(L,Tbar,'g','linewidth',2)
plot(L,Tmax,'k','color',[.5 0 .5],'linewidth',2)
plot(L,Tact,'r','linewidth',2)
legend([h1 h2 h3 h4],'Moreno data','Vanh pass','Vanh maxact','Vanh act','location','northwest')

V = (L / (2 * pi)).^2 * l * pi; 
Ppas = 2 * pi * Tpas ./ L; 
Pbar = 2 * pi * Tbar ./ L; 
Pmax = 2 * pi * Tmax ./ L; 
Pact = 2 * pi * Tact ./ L;

%% convert back to mmHg and mL 
V_mL = V * 1e6; 

Ppas_mmHg = Ppas / mmHg2SI; 
Pbar_mmHg = Pbar / mmHg2SI; 
Pmax_mmHg = Pmax / mmHg2SI; 
Pact_mmHg = Pact / mmHg2SI; 

Pd_Moreno_mmHg = Pd_Moreno / mmHg2SI; 
Vd_Moreno_mL = Vd_Moreno * 1e6; 

Pd_pass_Vanh_mmHg   = Pd_pass_Vanh / mmHg2SI; 
Pd_maxact_Vanh_mmHg = Pd_maxact_Vanh / mmHg2SI; 
Pd_act_Vanh_mmHg    = Pd_act_Vanh / mmHg2SI; 

Vd_pass_sca_mL   = Vd_pass_sca * 1e6; 
Vd_maxact_sca_mL = Vd_maxact_sca * 1e6; 

V_0_mL = V_0 * 1e6; 
V_M_mL = V_M * 1e6; 

figure(7) 
clf 
hold on 
plot([0 V_M_mL+1],zeros(2,1),'k')
plot([0 V_M_mL+1],P_n * ones(2,1),'k:')
plot([0 V_M_mL+1],P_M * ones(2,1),'k:')
plot(zeros(2,1),[-20 50],'k')
plot(V_0_mL * ones(2,1),[-20 50],'k:')
plot(V_n * ones(2,1),[-20 50],'k:')
plot(V_M_mL * ones(2,1),[-20 50],'k:')
l1 = plot(Vd_Moreno_mL,Pd_Moreno_mmHg,'mo','MarkerSize',5,'linewidth',2);
l2 = plot(Vd_pass_sca_mL,Pd_pass_Vanh_mmHg,'bo','MarkerSize',5,'linewidth',2);
l3 = plot(Vd_maxact_sca_mL,Pd_maxact_Vanh_mmHg,'ko','color',[.5 0 .5],'MarkerSize',5,'linewidth',2);
l4 = plot(Vd_pass_sca_mL,Pd_act_Vanh_mmHg,'ro','MarkerSize',5,'linewidth',2);
plot(V_mL,Ppas_mmHg,'b','linewidth',2)
plot(V_mL,Pbar_mmHg,'g','linewidth',2)
plot(V_mL,Pmax_mmHg,'k','color',[.5 0 .5],'linewidth',2)
plot(V_mL,Pact_mmHg,'r','linewidth',2)
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
set(gca,'FontSize',15)
legend([l1 l2 l3 l4],'Moreno data','Vanh pass','Vanh maxact','Vanh act','location','northwest')
xlim([0 V_M_mL+1])
ylim([-20 50])









