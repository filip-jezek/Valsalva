%% plots the baseline results
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
dl_n = dymload('../../Results/CardiovascularSystem.mat')

dl_avst = dymload('../../Results/imp_avSt.mat')
dl_avre = dymload('../../Results/imp_avRe.mat')

%%
color_b = [28, 108, 200]/255;
color_r = [238, 46, 47]/255;
color_g = [0, 140, 72]/255;
color_m = [226, 113, 199]/255;
color_lb = [182, 226, 255]/255;
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
mlPmin2SI = 1/1000/60;
%%
time_n = dymget(dl_n, 'Time');
d = 1.6;
t_interval_n = [28.2, 28.2 + d]; % interval in seconds
td = t_interval_n(2) - t_interval_n(1)
i_int_n = (time_n >= t_interval_n(1) & time_n <= t_interval_n(2));
t_n = time_n - t_interval_n(1);

pb = dymget(dl_n, 'brachial_pressure')/mmHg2SI;
pbm = dymget(dl_n, 'brachial_pressure_mean')/mmHg2SI;
hr = dymget(dl_n, 'HR')/bpm2SI;
plv = dymget(dl_n, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
vlv = dymget(dl_n, 'V_LV')/ml2SI;
pow_lv_n = dymget(dl_n, 'heartComponent.ventricles.power_LV');


time_avst = dymget(dl_avst, 'Time');
t_interval_avst = [298.28, 298.16 + d]; % interval in seconds
i_int_avst = (time_avst >= t_interval_avst(1) & time_avst <= t_interval_avst(2));
t_avst = time_avst - t_interval_avst(1);

pbm_avst = dymget(dl_avst, 'brachial_pressure_mean')/mmHg2SI;
pb_avst = dymget(dl_avst, 'brachial_pressure')/mmHg2SI;
hr_avst = dymget(dl_avst, 'HR')/bpm2SI;
plv_avst = dymget(dl_avst, 'P_LV')/mmHg2SI;
vlv_avst = dymget(dl_avst, 'V_LV')/ml2SI;
pascaor = dymget(dl_avst, 'heartComponent.sa.pressure')/mmHg2SI;
pow_lv_avst = dymget(dl_avst, 'heartComponent.ventricles.power_LV');
% figure(101);clf;hold on;
% plot(time_avst, pascaor)
% plot(time_avst, plv_avst)
disp(['Pressure drop across AV:' num2str(max(plv_avst(i_int_avst)) - max(pascaor(i_int_avst)) ) ])
disp(['power N:' num2str(max(pow_lv_n(i_int_n))) 'power AVST:' num2str(max(pow_lv_avst(i_int_avst)))])

time_avre = dymget(dl_avre, 'Time');
t_interval_avre = [298.35, 300]; % interval in seconds
i_int_avre = (time_avre >= t_interval_avre(1) & time_avre <= t_interval_avre(2));
t_avre = time_avre - t_interval_avre(1);

pbm_avre = dymget(dl_avre, 'brachial_pressure_mean')/mmHg2SI;
pb_avre = dymget(dl_avre, 'brachial_pressure')/mmHg2SI;
hr_avre = dymget(dl_avre, 'HR')/bpm2SI;
co_avre = dymget(dl_avre, 'CO')*60000;
plv_avre = dymget(dl_avre, 'P_LV')/mmHg2SI;
vlv_avre = dymget(dl_avre, 'heartComponent.ventricles.V_LV')/ml2SI;
pow_lv_avre = dymget(dl_avre, 'heartComponent.ventricles.power_LV');
regf = dymget(dl_avre, 'regurgitation_frac');
disp(['Regurg ' num2str(max(regf(i_int_avre)))] )
disp(['HR avst' num2str(mean(hr_avst(end)))])
disp(['HR avre' num2str(mean(hr_avre(end)))])

%%
fig1 = figure(1);clf;
fs = 8;
set(gcf, 'DefaultAxesFontSize', fs);

s_a1 = subplot(3, 2, 1);
hold on;
title('A: PV loops of AV stenosis ', 'FontSize', fs + 2);
plot(vlv(i_int_n), plv(i_int_n), 'Color',color_b, 'LineWidth', 1);
plot(vlv_avst(i_int_avst), plv_avst(i_int_avst), 'Color',color_r,  'LineWidth', 1);
xlim([20, 200])
ylim([0, 150])
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
% yyaxis right
% ylabel('W')
% ylim([0, 3.6])
% plot([50, 200], pow_lv_n(find(time_n > t_interval_n(1), 2)), 'b--', 'LineWidth', 2);
% plot([50, 200], pow_lv_avst(find(time_avst > t_interval_avst(1), 2)), 'r--', 'LineWidth', 2);
% leg = legend('LV N', 'LV AvSt', 'LV power N', 'LV power AvSt', 'Location', 'SouthWest')
leg = legend('LV N', 'LV AvSt', 'Location', 'SouthWest')
leg.ItemTokenSize = [10, 10];

s_a2 = subplot(3, 2, 2); 
hold on;
title('B: Arterial pressure of AV stenosis', 'FontSize', fs + 2);
plot(t_n(i_int_n), pb(i_int_n), 'Color',color_b, 'LineWidth', 1);
plot(t_n(i_int_n), pbm(i_int_n), '--','Color',color_b, 'LineWidth', 0.5);

plot(t_avst(i_int_avst), pb_avst(i_int_avst), 'Color',color_r,  'LineWidth', 1);
plot(t_avst(i_int_avst), pbm_avst(i_int_avst), '--', 'Color',color_r, 'LineWidth', 0.5);

% plot(t_n(i_int_n), hr(i_int_n), 'b:', 'LineWidth', 1.5);
% plot(t_avst(i_int_avst), hr_avst(i_int_avst), 'r:', 'LineWidth', 1.5);
% leg = legend('BP N', 'BPm N', 'BP AvSt', 'BPm AvSt', 'HR, N', 'HR AvSt', 'Location', 'SouthWest')
leg = legend('PA N', 'PA mean N', 'PA AvSt', 'PA mean AvSt', 'Location', 'best')
leg.ItemTokenSize = [10, 10];
xlabel('t (s)')
xlim([0, d]);
ylim([70, 120])
s_a2.Clipping = 'off';


s_a3 = subplot(3, 2, 3);
hold on;
title('C: PV loops of AV insufficiency', 'FontSize', fs + 2);
plot(vlv(i_int_n), plv(i_int_n), 'Color',color_b, 'LineWidth', 1);
plot(vlv_avre(i_int_avre), plv_avre(i_int_avre), 'Color',color_r,  'LineWidth', 1);
xlim([20, 200])
ylim([0, 150])
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
% yyaxis right
% ylabel('W')
% ylim([0, 3.6])
% plot([50, 250], pow_lv_n(find(time_n > t_interval_n(1), 2)), 'b--', 'LineWidth', 2);
% plot([50, 250], pow_lv_avre(find(time_avre > t_interval_avre(1), 2)), 'r--', 'LineWidth', 2);
% leg = legend('LV N', 'LV AvSt', 'LV power N', 'LV power AvSt', 'Location', 'NorthEast')
leg = legend('LV N', 'LV AvIns', 'Location', 'NorthEast')
leg.ItemTokenSize = [10, 10];

s_a4 = subplot(3, 2, 4); 
hold on;
title('D: Arterial pressure of AV insufficiency', 'FontSize', fs + 2);
plot(t_n(i_int_n), pb(i_int_n), 'Color',color_b, 'LineWidth', 1);
plot(t_n(i_int_n), pbm(i_int_n), '--', 'Color',color_b, 'LineWidth', 0.5);

plot(t_avre(i_int_avre), pb_avre(i_int_avre), 'Color',color_r,  'LineWidth', 1);
plot(t_avre(i_int_avre), pbm_avre(i_int_avre), '--', 'Color',color_r, 'LineWidth', 0.5);

% plot(t_n(i_int_n), hr(i_int_n), 'b:', 'LineWidth', 1.5);
% plot(t_avre(i_int_avre), hr_avre(i_int_avre), 'r:', 'LineWidth', 1.5);
leg = legend('PA N', 'PA mean N', 'PA AvIns', 'PA mean AvIns', 'HR N', 'HR AvIns', 'Location', 'best')
leg.ItemTokenSize = [10, 10];
xlabel('t (s)')
xlim([0, d]);
ylim([-10, 120])
yticks([0, 20, 40, 60, 80, 100, 120])
s_a4.Clipping = 'off';

%% Added art stiffening here

dl_n = dymload('../../Results/CardiovascularSystem.mat');
% dl_arst_m = '../../Results2/imp_arst.mat'
% system("c:\Program Files\Dymola 2021x\bin\dsres2sdf.exe ../../Results2/imp_arst.mat")
dl_arst = dymload('../../Results/imp_arst.mat')

time_n = dymget(dl_n, 'Time');
d = 1.6;
t_interval_n = [28.2, 28.2 + d]; % interval in seconds
td = t_interval_n(2) - t_interval_n(1)
i_int_n = (time_n >= t_interval_n(1) & time_n <= t_interval_n(2));
t_n = time_n - t_interval_n(1);

pb = dymget(dl_n, 'brachial_pressure')/mmHg2SI;
pbm = dymget(dl_n, 'brachial_pressure_mean')/mmHg2SI;
hr = dymget(dl_n, 'HR')/bpm2SI;
plv = dymget(dl_n, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
vlv = dymget(dl_n, 'V_LV')/ml2SI;
pow_lv_n = dymget(dl_n, 'heartComponent.ventricles.power_LV');


time_arst = dymget(dl_arst, 'Time');
t_interval_arst = [397.85, 398.16 + d]; % interval in seconds
i_int_arst = (time_arst >= t_interval_arst(1) & time_arst <= t_interval_arst(2));
t_arst = time_arst - t_interval_arst(1);

pbm_arst = dymget(dl_arst, 'brachial_pressure_mean')/mmHg2SI;
pb_arst = dymget(dl_arst, 'brachial_pressure')/mmHg2SI;
hr_arst = dymget(dl_arst, 'HR')/bpm2SI;
plv_arst = dymget(dl_arst, 'P_LV')/mmHg2SI;
vlv_arst = dymget(dl_arst, 'V_LV')/ml2SI;
pow_lv_arst = dymget(dl_arst, 'heartComponent.ventricles.power_LV');

disp(['PP arst' num2str(max(pb_arst(i_int_arst)) - min(pb_arst(i_int_arst)) )])
disp(['PAm diff' num2str(mean(pbm(i_int_n)) - mean(pbm_arst(i_int_arst)) )])
disp(['HR' num2str(mean(hr_arst(end)))])
%%
fig1 = figure(1);
fs = 8;
set(gcf, 'DefaultAxesFontSize', fs);

s_a1 = subplot(3, 2, 5);cla;hold on;
title('E: PV loops of Art. stiffening', 'FontSize', fs + 2);
% subplot(1, 2, 1)
plot(vlv(i_int_n), plv(i_int_n), 'Color',color_b, 'LineWidth', 1);
plot(vlv_arst(i_int_arst), plv_arst(i_int_arst), 'Color',color_r, 'LineWidth', 1);
% plot(t_arst(i_int_arst), plv_arst(i_int_arst), 'Color',color_b, 'LineWidth', 1);
leg = legend('LV N', 'LV ArSt', 'Location', 'NorthEast')
leg.ItemTokenSize = [10, 10];
xlim([20, 200])
ylim([0, 150])
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')

s_a1 = subplot(3, 2, 6);cla;hold on;
title('F: Arterial pressure of Art. stiffening', 'FontSize', fs + 2);
plot(t_n(i_int_n), pb(i_int_n), 'Color',color_b, 'LineWidth', 1);
plot(t_n(i_int_n), pbm(i_int_n), '--','Color',color_b, 'LineWidth', 0.5);

plot(t_arst(i_int_arst), pb_arst(i_int_arst), 'Color',color_r, 'LineWidth', 1);

plot(t_arst(i_int_arst), pbm_arst(i_int_arst), '--', 'Color',color_r,'LineWidth', 0.5);

% plot(t_n(i_int_n), hr(i_int_n), 'r:', 'LineWidth', 1);
% plot(t_arst(i_int_arst), hr_arst(i_int_arst), 'r', 'LineWidth', 0.5);


leg = legend('PA N', 'PA mean N', 'PA ArSt', 'PA mean ArSt', 'Location', 'SouthWest')
leg.ItemTokenSize = [10, 10];
ylabel('Pressure (mmHg)')
xlabel('t (s)')
xlim([0, d]);

%% save figure
tw = 17;
th = 17;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th], 'PaperUnits', 'Centimeters', 'PaperSize', [tw, th])

%%
% dsasdgfa
exportgraphics(gcf,'fig_R_ValvePath.png','Resolution',300)
exportgraphics(gcf,'fig_R_ValvePath.pdf', 'ContentType','vector')
% print(gcf,'fig_R_base_300.png', '-dpng', '-r300')
% saveas(gcf,'fig1.svg')