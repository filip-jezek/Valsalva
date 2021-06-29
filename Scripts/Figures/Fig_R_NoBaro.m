%% plots the baseline results
dl_vm = '../../Results/CVS_VMNoBaro.mat'
dl_ti = '../../Results/CVS_TiltNoBaro.mat'
% c:\Program Files\Dymola 2021x\bin\dsres2sdf.exe dl_n
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
dl_vm = dymload(dl_vm);
dl_ti = dymload(dl_ti);

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
cutofftime = 5;
decimateFactor = 10
time_vm = decimate(dymget(dl_vm, 'Time'), decimateFactor);
d = 1.6;
t_interval_n = [58.07, 58.2 + d]; % interval in seconds
td = t_interval_n(2) - t_interval_n(1)
i_int_n = (time_vm >= t_interval_n(1) & time_vm <= t_interval_n(2));
t_n = time_vm - t_interval_n(1);

time_ti = decimate(dymget(dl_vm, 'Time'), decimateFactor, 10);
pb_vm = decimate(dymget(dl_vm, 'brachial_pressure')/mmHg2SI, decimateFactor, 10); 
pbs_vm = decimate(dymget(dl_vm, 'brachial_pressure_systolic')/mmHg2SI, decimateFactor, 10);
pbs_vm(time_vm < cutofftime) = pbs_vm(find(time_ti >= cutofftime, 1));
pbd_vm = decimate(dymget(dl_vm, 'brachial_pressure_diastolic')/mmHg2SI, decimateFactor, 10);
pbd_vm(time_vm < cutofftime) = pbd_vm(find(time_ti >= cutofftime, 1));
pbm_vm = decimate(dymget(dl_vm, 'brachial_pressure_mean')/mmHg2SI, decimateFactor, 10);
pbm_vm(time_vm < cutofftime) = pbm_vm(find(time_ti >= cutofftime, 1));
co_vm = decimate(dymget(dl_vm, 'CO')/ml2SI*60/1000, decimateFactor, 10);
co_vm(time_vm < cutofftime) = co_vm(find(time_ti >= cutofftime, 1));


time_ti = decimate(dymget(dl_ti, 'Time'), decimateFactor, 10);
pb_ti = decimate(dymget(dl_ti, 'brachial_pressure')/mmHg2SI, decimateFactor, 10);
pbs_ti = decimate(dymget(dl_ti, 'brachial_pressure_systolic')/mmHg2SI, decimateFactor, 10);
pbs_ti(time_ti < cutofftime) = pbs_ti(find(time_ti >= cutofftime, 1));
pbd_ti = decimate(dymget(dl_ti, 'brachial_pressure_diastolic')/mmHg2SI, decimateFactor, 10);
pbd_ti(time_ti < cutofftime) = pbd_ti(find(time_ti >= cutofftime, 1));
pbm_ti = decimate(dymget(dl_ti, 'brachial_pressure_mean')/mmHg2SI, decimateFactor, 10);
pbm_ti(time_ti < cutofftime) = pbm_ti(find(time_ti >= cutofftime, 1));
co_ti = decimate(dymget(dl_ti, 'CO')/ml2SI*60/1000, decimateFactor, 10);
co_ti(time_ti < cutofftime) = co_ti(find(time_ti >= cutofftime, 1));

%%
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 8, 'defaultLineLineWidth',1.0);

s_a1 = subplot(2, 1, 1);
hold on;
title('A: VM with impaired baroreflex');
% title('\fontsize{16}A: \fontsize{12}VM with impaired baroreflex');
% plot(time_vm, pb_vm, 'b', 'LineWidth', 0.5)
fill([time_vm; flipud(time_vm)], [pbs_vm; flipud(pbd_vm)], color_lb,'EdgeColor',[182/255, 226/255, 1])
plot(time_vm, pbm_vm, 'Color', color_b,'LineWidth', 1)
xlim([0,60])
ylim([20, 140])
xlabel('t (s)')
ylabel('Pressure (mmHg)')

yyaxis right;
plot(time_vm, co_vm,'--','Color', color_m, 'LineWidth', 1.5)
ylim([0 8])
s_a1.YAxis(2).Color = [0.4940, 0.1840, 0.5560];
% set(gca,'ytick',[])
ylabel('CO (L/min)')
leg = legend('PA (mmHg)', 'PA mean', 'CO (L/min)', 'Location', 'SouthWest');
leg.ItemTokenSize = [10, 2];
g = gca; g.YAxis(2).Color = color_m;

s_a2 = subplot(2, 1, 2);
hold on;
% title('\fontsize{16}B: \fontsize{12}Tilt with impaired baroreflex');
title('B: HUT with impaired baroreflex');
% plot(time_ti, pb_ti, 'b', 'LineWidth', 0.5)

fill([time_ti; flipud(time_ti)], [pbs_ti; flipud(pbd_ti)], color_lb,'EdgeColor',[182/255, 226/255, 1])
plot(time_ti, pbm_ti, 'Color', color_b,'LineWidth', 1.5)
xlim([min(time_ti), max(time_ti)])
ylim([60, 120])
xlabel('t (s)')
% set(gca,'ytick',[])
ylabel('Pressure (mmHg)')

yyaxis right;
plot(time_ti, co_ti,'--' ,'Color', color_m, 'LineWidth', 1.5)
ylim([0 8])
s_a2.YAxis(2).Color = color_m;
ylabel('CO (L/min)')

leg = legend('PA (mmHg)', 'PA mean', 'CO (L/min)', 'Location', 'SouthEast');
leg.ItemTokenSize = [15, 2];

%% save figure
tw = 7;
th = 12;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th], 'PaperUnits', 'Centimeters', 'PaperSize', [tw, th])
drawnow()
%%
pos = get(s_a2, 'Position');

w = 0.75
h = 0.4
% In row alignment
% set(s_a1, 'Position', [0.08, 0.17, w, h])
% set(s_a2, 'Position', [0.5, 0.17, w, h])
% column alignment


% set(s_a1, 'Position', [0.08, 0.17 + h, w, h])
% set(s_a2, 'Position', [0.08, 0.17, w, h])

% s_a1.Title.Units = 'normalized';
% s_a1.Title.HorizontalAlignment = 'left';
% s_a1.Title.VerticalAlignment = 'bottom';
% s_a1.Title.Position = [0, 1.02, 0];
% 
% s_a2.Title.Units = 'normalized';
% s_a2.Title.HorizontalAlignment = 'left';
% s_a2.Title.VerticalAlignment = 'bottom';
% s_a2.Title.Position = [0, 1.02, 0];
%%
exportgraphics(gcf,'fig_R_NoBaro.png','Resolution',300)
exportgraphics(gcf,'fig_R_NoBaro.pdf', 'ContentType','vector')
% print(gcf,'fig_R_base_300.png', '-dpng', '-r300')
% saveas(gcf,'fig1.svg')