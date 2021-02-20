%% plots the baseline results
dl_n = '../../Results2/CardiovascularSystem.mat'
dl_avst = '../../Results2/imp_avSt.sdf'
dl_avre = '../../Results2/imp_avRe.sdf'
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
dl_n = dymload(dl_n)
%%
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
%%
time_n = dymget(dl_n, 'Time');
d = 1.6;
t_interval_n = [58.2, 58.2 + d]; % interval in seconds
td = t_interval_n(2) - t_interval_n(1)
i_int_n = (time_n >= t_interval_n(1) & time_n <= t_interval_n(2));
t_n = time_n - t_interval_n(1);

pb = dymget(dl_n, 'brachial_pressure')/mmHg2SI;
pbm = dymget(dl_n, 'brachial_pressure_mean')/mmHg2SI;
hr = dymget(dl_n, 'HR')/bpm2SI;
plv = dymget(dl_n, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
vlv = dymget(dl_n, 'V_LV')/ml2SI;
pow_lv_n = dymget(dl_n, 'heartComponent.ventricles.power_LV');


time_avst = h5read(dl_avst, '/Time');
t_interval_avst = [298.16, 298.16 + d]; % interval in seconds
i_int_avst = (time_avst >= t_interval_avst(1) & time_avst <= t_interval_avst(2));
t_avst = time_avst - t_interval_avst(1);

pbm_avst = h5read(dl_avst, '/brachial_pressure_mean')/mmHg2SI;
pb_avst = h5read(dl_avst, '/brachial_pressure')/mmHg2SI;
hr_avst = h5read(dl_avst, '/HR')/bpm2SI;
plv_avst = h5read(dl_avst, '/P_LV')/mmHg2SI;
vlv_avst = h5read(dl_avst, '/heartComponent/ventricles/V_LV')/ml2SI;
pow_lv_avst = h5read(dl_avst, '/heartComponent/ventricles/power_LV');

time_avre = h5read(dl_avre, '/Time');
t_interval_avre = [298.16, 298.16 + d]; % interval in seconds
i_int_avre = (time_avre >= t_interval_avre(1) & time_avre <= t_interval_avre(2));
t_avre = time_avre - t_interval_avre(1);

pbm_avre = h5read(dl_avre, '/brachial_pressure_mean')/mmHg2SI;
pb_avre = h5read(dl_avre, '/brachial_pressure')/mmHg2SI;
hr_avre = h5read(dl_avre, '/HR')/bpm2SI;
co_avre = h5read(dl_avre, '/CO')*60000;
plv_avre = h5read(dl_avre, '/P_LV')/mmHg2SI;
vlv_avre = h5read(dl_avre, '/heartComponent/ventricles/V_LV')/ml2SI;
pow_lv_avre = h5read(dl_avre, '/heartComponent/ventricles/power_LV');

%%
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 10);

s_a1 = subplot(2, 2, 1);
hold on;
title('A: PV loop of AV stenosis ');
plot(vlv(i_int_n), plv(i_int_n), 'b', 'LineWidth', 1);
plot(vlv_avst(i_int_avst), plv_avst(i_int_avst), 'r', 'LineWidth', 1);
xlim([50, 180])
ylim([0, 180])
ylabel('Pressure (mmHg)')
xlabel('Volume (ml)')
% yyaxis right
% ylabel('W')
% ylim([0, 3.6])
% plot([50, 200], pow_lv_n(find(time_n > t_interval_n(1), 2)), 'b--', 'LineWidth', 2);
% plot([50, 200], pow_lv_avst(find(time_avst > t_interval_avst(1), 2)), 'r--', 'LineWidth', 2);
% leg = legend('LV N', 'LV AvSt', 'LV power N', 'LV power AvSt', 'Location', 'SouthWest')
leg = legend('LV N', 'LV AvSt', 'Location', 'SouthWest')
leg.ItemTokenSize = [10, 10];

s_a2 = subplot(2, 2, 2); 
hold on;
title('B: BP and HR of AV stenosis');
plot(t_n(i_int_n), pb(i_int_n), 'b', 'LineWidth', 1);
plot(t_n(i_int_n), pbm(i_int_n), 'b--', 'LineWidth', 0.5);

plot(t_avst(i_int_avst), pb_avst(i_int_avst), 'r', 'LineWidth', 1);
plot(t_avst(i_int_avst), pbm_avst(i_int_avst), 'r--', 'LineWidth', 0.5);

% plot(t_n(i_int_n), hr(i_int_n), 'b:', 'LineWidth', 1.5);
% plot(t_avst(i_int_avst), hr_avst(i_int_avst), 'r:', 'LineWidth', 1.5);
% leg = legend('BP N', 'BPm N', 'BP AvSt', 'BPm AvSt', 'HR, N', 'HR AvSt', 'Location', 'SouthWest')
leg = legend('BP N', 'BPm N', 'BP AvSt', 'BPm AvSt', 'Location', 'SouthWest')
leg.ItemTokenSize = [10, 10];
xlabel('t (s)')
xlim([0, d]);
ylim([50, 120])
s_a2.Clipping = 'off';


s_a3 = subplot(2, 2, 3);
hold on;
title('C: PV loop of AV insufficiency');
plot(vlv(i_int_n), plv(i_int_n), 'b', 'LineWidth', 1);
plot(vlv_avre(i_int_avre), plv_avre(i_int_avre), 'r', 'LineWidth', 1);
xlim([50, 250])
ylim([0, 120])
ylabel('Pressure (mmHg)')
xlabel('Volume (ml)')
% yyaxis right
% ylabel('W')
% ylim([0, 3.6])
% plot([50, 250], pow_lv_n(find(time_n > t_interval_n(1), 2)), 'b--', 'LineWidth', 2);
% plot([50, 250], pow_lv_avre(find(time_avre > t_interval_avre(1), 2)), 'r--', 'LineWidth', 2);
% leg = legend('LV N', 'LV AvSt', 'LV power N', 'LV power AvSt', 'Location', 'NorthEast')
leg = legend('LV N', 'LV AvSt', 'Location', 'NorthEast')
leg.ItemTokenSize = [10, 10];

s_a4 = subplot(2, 2, 4); 
hold on;
title('D: BP and HR of AV insufficiency');
plot(t_n(i_int_n), pb(i_int_n), 'b', 'LineWidth', 1);
plot(t_n(i_int_n), pbm(i_int_n), 'b--', 'LineWidth', 0.5);

plot(t_avre(i_int_avre), pb_avre(i_int_avre), 'r', 'LineWidth', 1);
plot(t_avre(i_int_avre), pbm_avre(i_int_avre), 'r--', 'LineWidth', 0.5);

% plot(t_n(i_int_n), hr(i_int_n), 'b:', 'LineWidth', 1.5);
% plot(t_avre(i_int_avre), hr_avre(i_int_avre), 'r:', 'LineWidth', 1.5);
leg = legend('BP N', 'BPm N', 'BP AvSt', 'BPm AvSt', 'HR, N', 'HR AvSt', 'Location', 'SouthWest')
leg.ItemTokenSize = [10, 10];
xlabel('t (s)')
xlim([0, d]);
ylim([0, 120])
s_a4.Clipping = 'off';

%% save figure
tw = 17;
th = 14;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th], 'PaperUnits', 'Centimeters', 'PaperSize', [tw, th])
%%
exportgraphics(gcf,'fig_R_ValvePath.png','Resolution',300)
exportgraphics(gcf,'fig_R_ValvePath.pdf', 'ContentType','vector')
% print(gcf,'fig_R_base_300.png', '-dpng', '-r300')
% saveas(gcf,'fig1.svg')