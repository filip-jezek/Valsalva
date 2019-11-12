
%% Load data

load Moreno70_data.mat
load Vanhoutte69_data.mat

%% Formulation 

l      = .02; %m 
phibar = 0.2; 

%Zero pressure/tension values at phi = 0 (Combining the Vanhoutte and
%Moreno passive data)
V0 = max(Vd_pass_Vanh)/3;       %mL
L0 = 2*pi*sqrt(V0*1e-6/l/pi);   %m

%Nominal point at phi = 0.2
Pn = 10*133;                    %Pa
Vn = 2*V0;                      %mL
Ln = 2*pi*sqrt(Vn*1e-6/l/pi);   %m
Tn = Pn*Ln/(2*pi);              %N

%Maximal pressure and volume point at phi = 1
PM = 40*133;                    %Pa
VM = 3*V0 + V0;                 %mL
LM = 2*pi*sqrt(VM*1e-6/l/pi);   %m
TM = PM*LM/(2*pi);              %N

%Functional components
f = @(L) L.*(L - L0)/L0;                        %m
g = @(L) L0*(exp( 11.5*(L - L0)/L0 ) - 1);      %m
h = @(L) (L - L0).*exp( -1.5*((L - L0)/L0).^2); %m

a = 1.2*133;                    %Pa
b = 6e-5*133;                   %Pa
c = 26.25*133;                  %Pa

%Tensions
T_P      = @(L) a*f(L) + b*g(L);
T_A      = @(L) c*h(L); 
T_nom    = @(L) T_P(L) + phibar*T_A(L);
T_maxact = @(L) T_P(L) + T_A(L);

%Circumferential length vector (m)
L = 1e-6:1e-6:LM; 

%Volume (mL)
V = pi.*(L/(2*pi)).^2*l/1e-6;
%Relative volume (dimensionless)
RelV = (V - V0)/V0;

%Pressures
P_nom    = @(L) T_nom(L)*(2*pi)./L; 
P_pas    = @(L) T_P(L)*(2*pi)./L; 
P_act    = @(L) T_A(L)*(2*pi)./L; 
P_maxact = @(L) T_maxact(L).*(2*pi)./L; 

%% Plots 

minT = min([min(T_nom(L)), min(T_P(L)), min(T_A(L))]); 
maxT = max([max(T_nom(L)), max(T_P(L)), max(T_A(L))]); 

%{
Color map

Curves: 
Nominal             Green
Passive             Blue
Maximally active    Purple  [.5 0 .5]
Active only         Red 

Points: 
Nominal             Gray    [.5 .5 .5]
Zero pressure       Black
Maximal pressure    Magenta
Max act shift       Orange  [1 .65 0]
%}

figure(1)
clf
hold on 
%axes
plot([min(L) max(L)], zeros(2,1), 'k')
%points
plot(L0*ones(2,1), [minT TM],'k:')
plot(Ln*ones(2,1), [minT TM],':','Color',[.5 .5 .5])
plot([min(L) max(L)], Tn*ones(2,1),':','Color',[.5 .5 .5])
plot(LM*ones(2,1), [minT TM],'m:')
plot([min(L) max(L)], TM*ones(2,1),'m:')
%curves
T3 = plot(L, T_maxact(L),'Color',[.5 0 .5],'linewidth',2);
T2 = plot(L, T_P(L),'b','linewidth',2);
T4 = plot(L, T_A(L),'r','linewidth',2);
T1 = plot(L, T_nom(L),'g','linewidth',2); 
xlabel('Circumferential length (m)')
xlim([min(L) max(L)])
ylim([minT-1 TM+1])
ylabel('Tension (N)') 
legend([T1 T2 T3 T4],... 
    '\phi = 0.2', ...
    '\phi = 0',...
    '\phi = 1',...
    'Active alone',... 
    'location','northwest') 

print -dpng TvsL.png

%V vs P curve
figure(2) 
clf
hold on
%axes
plot([-20 50], zeros(2,1),'k')
plot(zeros(2,1),[0 VM*1.1],'k')
%points
plot([-20 50], V0*ones(2,1),'k:')
plot([-20 50], Vn*ones(2,1),':','Color',[.5 .5 .5])
plot([-20 50], VM*ones(2,1),'m:')
plot(Pn*ones(2,1), [0 VM*1.1], ':', 'Color',[.5 .5 .5])
plot(PM*ones(2,1), [0 VM*1.1], 'm:')
%curve
V1 = plot(P_nom(L)/133, V, 'g','linewidth',2);
V2 = plot(P_pas(L)/133, V, 'b','linewidth',2);
V3 = plot(P_maxact(L)/133, V, 'Color',[.5 0 .5], 'linewidth',2);
V4 = plot(P_act(L)/133, V, 'r', 'linewidth',2); 
xlabel('Pressure (mmHg)')
ylabel('Volume (mL)')
ylim([0 VM*1.1])
xlim([-20 50])
legend([V1 V2 V3 V4],'\phi = 0.2','\phi = 0', '\phi = 1','Active alone','location','northwest')

print -dpng VvsP.png


Pd_actonly_Vanh = Pd_maxact_Vanh - Pd_pass_Vanh; 
Vd_actonly_Vanh = Vd_maxact_Vanh + V0; 
Ld_actonly_Vanh = 2*pi*sqrt(Vd_actonly_Vanh/pi/l); 
Td_actonly_Vanh = Pd_actonly_Vanh.*Ld_actonly_Vanh/(2*pi); 

%P vs V curve
figure(3) 
clf
hold on
%axes
plot(zeros(2,1),[-20 50], 'k')
plot([0 VM*1.1],zeros(2,1),'k')
%points
plot(V0*ones(2,1),[-20 50], 'k:')
plot(Vn*ones(2,1),[-20 50], ':','Color',[.5 .5 .5])
plot(VM*ones(2,1),[-20 50], 'm:')
plot([0 VM*1.1], Pn/133*ones(2,1),  ':', 'Color',[.5 .5 .5])
plot([0 VM*1.1], PM/133*ones(2,1),  'm:')
%data - with scaled volume
plot((V0 + Vd_pass_Vanh), Pd_pass_Vanh , 'o','Color',[0 0 .8],'linewidth',2,'Markersize',5)
plot((V0 + Vd_maxact_Vanh) , Pd_maxact_Vanh ,  'o','Color',[0.5 0 0.5],'linewidth',2,'Markersize',5)
plot(Vd_actonly_Vanh, Pd_actonly_Vanh,'o','Color',[.8 0 0],'linewidth',2,'Markersize',5)
%curves
P1 = plot(V, P_nom(L)/133,'g','linewidth',2);
P2 = plot(V, P_pas(L)/133, 'b', 'linewidth',2);
P3 = plot(V,P_maxact(L)/133,  'Color', [.5 0 .5], 'linewidth',2);
P4 = plot(V, P_act(L)/133, 'r', 'linewidth',2); 
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
ylim([-20 50])
xlim([0 VM*1.1])
legend([P1 P2 P3 P4],'\phi = 0.2','\phi = 0', '\phi = 1','Active alone','location','northwest')

print -dpng PvsV.png

%Rel V vs P curve
figure(4)
clf
hold on 
%axes
plot([-20 50], zeros(2,1),'k')
plot(zeros(2,1),[-1 3.5],'k')
%points
plot(Pn/133*ones(2,1),[-1 3.5],':','Color',[.5 .5 .5])
plot([-20 50], (Vn - V0)/V0*ones(2,1),':','Color',[.5 .5 .5])
plot(PM/133*ones(2,1),[-1 3.5],'m:')
plot([-20 50], (VM - V0)/V0*ones(2,1),'m:')
%data
plot(Pd_Moreno, Vd_sca_Moreno, 'o','Color',[.5 .5 .5],'linewidth',2,'Markersize',5)
%curves
r1 = plot(P_nom(L)/133, RelV, 'g','linewidth',2); 
r2 = plot(P_pas(L)/133, RelV, 'b','linewidth',2);
r3 = plot(P_maxact(L)/133, RelV, 'Color',[.5 0 .5],'linewidth',2);
r4 = plot(P_act(L)/133, RelV, 'r', 'linewidth',2);
xlabel('Pressure (mmHg)')
ylabel('Relative volume (V - V_0)/V_0')
xlim([-20 50])
ylim([-1 3.5])
%legend([r1 r2 r3],'\phi = 0.25','\phi = 0', '\phi = 1','location','northwest')
legend([r1 r2 r3 r4],'\phi = 0.2','\phi = 0', '\phi = 1','Active alone','location','northwest')

print -dpng RelVvsP.png







