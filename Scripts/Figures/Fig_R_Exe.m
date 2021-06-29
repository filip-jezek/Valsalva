%% plots the Exercise results
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
% system("c:\Program Files\Dymola 2021x\bin\dsres2sdf.exe imp_stepEx_normal")
datafile = '../../Results/imp_stepEx_normal.mat';

%% color defiinitons
color_b = [28, 108, 200]/255;
color_r = [238, 46, 47]/255;
color_g = [0, 140, 72]/255;
color_m = [226, 113, 199]/255;
color_lb = [182, 226, 255]/255;
color_s = [0.8, 0.8, 0.8];
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
%%

dl = dymload(datafile)
time = dymget(dl, 'Time');
t_interval = [0, 60]; % interval in seconds
td = t_interval(2) - t_interval(1);
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int = time >= 2; % get rid of the initial zeros

safezone = 400;
% %%
% g = hi.Groups(11)
% for i = 1:size(g.Groups, 1)
%     disp(num2str(i) + ":" + g.Groups(i).Name)
% end
%% EXERCISE VALUES
t = time;
pb = dymget(dl, 'brachial_pressure')/mmHg2SI;

pbs = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pbs(1:safezone) = pbs(safezone);
pbd = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pbd(1:safezone) = pbd(safezone);
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
pbm(1:safezone) = pbm(safezone);
psa = dymget(dl, 'heartComponent.sa.pressure')/mmHg2SI;

hr = dymget(dl, 'HR')/bpm2SI;
lsv = dymget(dl, 'SV')/ml2SI;
% rsv = dymget(dl, 'heartComponent/pulmonaryValve/SV')/ml2SI;
plv = dymget(dl, 'P_LV')/mmHg2SI;
prv = dymget(dl, 'heartComponent.ventricles.P_RV')/mmHg2SI;
vlv = dymget(dl, 'heartComponent.ventricles.V_LV')/ml2SI;
vrv = dymget(dl, 'heartComponent.ventricles.V_RV')/ml2SI;
vla = dymget(dl, 'V_la')/ml2SI;

pla = dymget(dl, 'P_LA')/mmHg2SI;

psv = dymget(dl, 'P_sv')/mmHg2SI;
ppv = dymget(dl, 'P_pv')/mmHg2SI;
ppa = dymget(dl, 'heartComponent.pa.pressure')/mmHg2SI;

pra = dymget(dl, 'heartComponent.tricuspidValve.q_in.pressure')/mmHg2SI;
vra = dymget(dl, 'heartComponent.ra.volume')/ml2SI;


tp = dymget(dl, 'thoracic_pressure')/mmHg2SI;
pdv = dymget(dl, 'SystemicComponent.femoral_vein_R34.port_a.pressure')/mmHg2SI;
psv = dymget(dl, 'P_sv')/mmHg2SI - tp(1);
ppv = dymget(dl, 'P_pv')/mmHg2SI - tp(1);
el = dymget(dl, 'Exercise.y_')*100;
qmv = dymget(dl, 'q_mv')/ml2SI/1000*60;

pow_lv = dymget(dl, 'heartComponent.ventricles.power_LV');
pow_rv = dymget(dl, 'heartComponent.ventricles.power_RV');

co = dymget(dl, 'CO')/ml2SI/1000*60;;
q_ex = dymget(dl, 'q_exercised_avg')/ml2SI/1000*60;

%% READ NORMAL FOR COMPARISON
datafile_n = '../../Results/CardiovascularSystem.mat'
dl_n = dymload(datafile_n)
time_n = dymget(dl_n, 'Time');
t_interval_n = [28.13, 28.13 + 1.6]; % interval in seconds
td_n = t_interval_n(2) - t_interval_n(1)
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int_n = (time_n >= t_interval_n(1) & time_n <= t_interval_n(2));

t_n = time_n - t_interval_n(1);
pb_n = dymget(dl_n, 'brachial_pressure')/mmHg2SI;
psa_n = dymget(dl_n, 'heartComponent.sa.pressure')/mmHg2SI;
plv_n = dymget(dl_n, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
pla_n = dymget(dl_n, 'heartComponent.mitralValve.q_in.pressure')/mmHg2SI;

psv_n = dymget(dl_n, 'P_sv')/mmHg2SI;
ppv_n = dymget(dl_n, 'P_pv')/mmHg2SI;


vlv_n = dymget(dl_n, 'V_LV')/ml2SI;
vla_n = dymget(dl_n, 'V_la')/ml2SI;

ppa_n = dymget(dl_n, 'heartComponent.pa.pressure')/mmHg2SI;
prv_n = dymget(dl_n, 'heartComponent.tricuspidValve.q_out.pressure')/mmHg2SI;
pra_n = dymget(dl_n, 'heartComponent.tricuspidValve.q_in.pressure')/mmHg2SI;

vrv_n = dymget(dl_n, 'heartComponent.ventricles.V_RV')/ml2SI;
vra_n = dymget(dl_n, 'heartComponent.ra.volume')/ml2SI;

qmv_n = dymget(dl_n, 'q_mv')/ml2SI/1000*60;


%% plot pv loops
ts = 35
te = 36.6
td = te - ts;
diffs = [1; find( diff(el) > 0)]
figure(3); clf;hold on;
ef = [];
edv = [];
pow_lvs = [];
pow_rvs = [];
co_s = [];
q_ex_s = [];

for i = 1:size(diffs)
%      i = 9

    i_s = find(time > time(diffs(i)) + ts, 1);
    i_e = find(time > time(diffs(i)) + te, 1);
%     t = time - time(i_s);
    i_int = [i_s:i_e];
    if i == 7 + 1
        t = time - time(i_s);
        i_int90 = i_int;
    end
%     plot(vlv(i_s:i_e), plv(i_s:i_e))
%     plot(vrv(i_s:i_e), prv(i_s:i_e))
%     plot(vla(i_s:i_e), pla(i_s:i_e))
    sv =  max(vlv(i_s:i_e)) -  min(vlv(i_s:i_e));
    ef = [ef, sv / max(vlv(i_s:i_e))];
    edv = [edv, min(vlv(i_s:i_e))];
    pow_lvs = [pow_lvs , pow_lv(i_s)];
    pow_rvs = [pow_rvs , pow_rv(i_s)];
    co_s = [co_s, mean(co(i_int))];
    q_ex_s = [q_ex_s, mean(q_ex(i_int))];
end

% plot([0:10:100], edv, 'b*-')
% plot([0:10:100], pow_rvs, 'r*-')
% plot([0:10:100], pow_lvs, 'b*-')
% plot([0:10:100], ef*10, 'b*-')
plot([0:10:100], co_s, 'b*-')
plot([0:10:100], q_ex_s, 'r*-')
% title('Ejection fraction on exercise level')

%% DRAW
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 8);

s_a1 = subplot(4, 2, 1);
hold on;
title('A: Systemic and LH pressures and volumes at 70% exercise');
plot(t(i_int90), plv(i_int90), 'Color', color_b, 'LineWidth', 1);
plot(t(i_int90), psa(i_int90), 'Color', color_r, 'LineWidth', 1);
plot(t(i_int90), pb(i_int90), 'Color', [0 0 0], 'LineWidth', 1);
plot(t(i_int90), pla(i_int90), 'Color', color_m, 'LineWidth', 1);
set(gca,'xtick',[])
leg = legend('P LV', 'P Asc Aor', 'PA (brach art)', 'P LA', 'Location', 'SouthEast')
leg.ItemTokenSize = [10, 150];
% legend('boxoff')
ylim([0, 220]);
ylabel('Pressure (mmHg)')
s_a1.Clipping = 'off';
xlim([0 td])
    
% volumes
s_a2 = subplot(4, 2, 3);hold on;
plot(t(i_int90), vlv(i_int90), 'Color', color_b, 'LineWidth', 1);
plot(t(i_int90), vla(i_int90), 'Color', color_m, 'LineWidth', 1);
leg = legend('V LV', 'V LA');
leg.ItemTokenSize = [10, 2];
set(gca,'xtickMode', 'auto')
ylim([0, 200]);
xlim([0 td])
xlabel('t (s)');
ylabel('Volume (mL)')
s_a1.Clipping = 'off';
% pos = get(s1, 'Position');
% set(s1, 'Position', [0.05, 0.5, 0.9,0.4])
% set(s2, 'Position', [0.05, 0.05, 0.9,0.4])

% B - venous plots
s_b1 = subplot(4, 2, 2);cla;
hold on;
title('B: Pulmonary and venous pressures');
plot(t_n(i_int_n), psv_n(i_int_n), 'b:','Color', color_b, 'LineWidth', 1);
plot(t_n(i_int_n), ppv_n(i_int_n), ':','Color', color_r, 'LineWidth', 1);

plot(t(i_int90), psv(i_int90), 'Color', color_b, 'LineWidth', 1);
plot(t(i_int90), ppv(i_int90), 'Color', color_r, 'LineWidth', 1);
s_b1.Clipping = 'off';
% set(gca,'xtick',[], 'ytick', [])
set(gca,'xtick',[], 'ytick', [5, 10, 15])
leg = legend('P VC N', 'P PV N', 'P VC 70% E', 'P PV 70% E')
leg.ItemTokenSize = [10, 2];
ylim([0, 20])
ylabel('Pressure (mmHg)')
xlim([0 td])

% Pulmonary arterial pressures
s_b2 = subplot(4, 2, 4);hold on;
plot(t_n(i_int_n), ppa_n(i_int_n), ':','Color', color_b, 'LineWidth', 1);

plot(t(i_int90), ppa(i_int90), 'Color', color_b, 'LineWidth', 1);
% set(gca,'xtickMode', 'auto', 'ytick', [])
ylim([0, 40]);
ylabel('Pressure (mmHg)')
xlim([0 td])
leg = legend('PPA N', 'PPA 70% E')
leg.ItemTokenSize = [10, 2];
xlabel('t (s)');

% pos = get(s1, 'Position');
% set(s1, 'Position', [0.05, 0.5, 0.9,0.4])
% set(s2, 'Position', [0.05, 0.05, 0.9,0.4])

% C
s_c = subplot(2, 2, 3);cla;hold on;
title('C: Ventricle power during step-up exercise')
plot([0:10:100], pow_lvs, '*-', 'Color', color_b,'LineWidth', 1)
plot([0:10:100], pow_rvs, '*-', 'Color', color_r, 'LineWidth', 1)
legend('LV', 'RV', 'Location', 'NorthWest')
xlabel('Exercise (% of max)');
ylabel('Power (W)');

% D
s_d = subplot(2, 2, 4);cla;hold on;
title('D: Cardiac output during step-up exercise')
% plot(vra(i_int), pra(i_int));
% plot(vla(i_int), pla(i_int));
plot([0:10:100], co_s, '*-', 'Color', color_m, 'LineWidth', 1)
plot([0:10:100], q_ex_s, '*--', 'Color', color_m, 'LineWidth', 1)
xlabel('Exercise (% of max)');
ylabel('Blood flow (L/min)');
leg = legend('Cardiac output', 'Flow through exercising tissues', 'Location', 'SouthEast')
leg.ItemTokenSize = [20, 150];
drawnow()

%% align the axis
drawnow()
pos = get(s_a2, 'Position');
tw = 17;
th = 14;
set(fig1, 'Units', 'Centimeters', 'Position', [2, 4, tw, th])

hm = 0.075;
vm = 0.08;
set(s_a1, 'Units', 'normalized', 'Position', [hm, 0.75, 0.5-1.2*hm,0.25 - vm/2])
set(s_a2, 'Units', 'normalized', 'Position', [hm, 0.5, 0.5-1.2*hm,0.25 - vm/2])

set(s_b1, 'Units', 'normalized', 'Position', [0.5 + hm, 0.75, 0.5-1.2*hm,0.25 - vm/2])
set(s_b2, 'Units', 'normalized', 'Position', [0.5 + hm, 0.5, 0.5-1.2*hm,0.25 - vm/2])

set(s_c, 'Units', 'normalized', 'Position', [hm, 1*vm, 0.5-1.2*hm,0.5 - 2.5*vm])
set(s_d, 'Units', 'normalized', 'Position', [0.5 + hm, 1*vm, 0.5-1.2*hm,0.5 - 2.5*vm])

% %% Draw boxes
% annotation(fig1,'rectangle', [0 0.5 0.5 0.5]);
% annotation(fig1,'rectangle', [0.5 0.5 0.5 0.5]);


%% save figure
exportgraphics(gcf,'fig_R_Exer.png','Resolution',200)
exportgraphics(gcf,'fig_R_Exer.pdf', 'ContentType','vector')
% tw = 17;
% th = 10;
% set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th], 'PaperUnits', 'Centimeters', 'PaperSize', [tw, th])
% print(gcf,'fig_R_Exer_200.png', '-dpng', '-r200')

% saveas(gcf,'fig1.svg')