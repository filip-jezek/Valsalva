%% plots the Tilt fresults
% import the dymload util
addpath('c:\Program Files\Dymola 2021x\Mfiles\dymtools\')

datafile = '../../Results/CVS_tiltable.mat';
dl = dymload(datafile)
%%
%% color defiinitons
color_b = [28, 108, 200]/255;
color_r = [238, 46, 47]/255;
color_g = [0, 140, 72]/255;
color_m = [226, 113, 199]/255;
color_lb = [182, 226, 255]/255;
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
%%
T = readtable('..\..\data\tilt\Wieling_dataset_PMID9640339.csv');
i_past = find(T.PAs_t > -10, 1)
i_padt = find(T.PAd_t > -10, 1)
i_hrt = find(T.HR_t > -10, 1)
i_coft = find(T.COf_t > -10, 1)
%%

time = dymget(dl, 'Time');
t_interval = [0, 130]; % interval in seconds
td = t_interval(2) - t_interval(1);
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int = [find(time >= t_interval(1), 1) : find(time >= t_interval(2), 1)] ; % get rid of the initial zeros
safezone = 400;


time = time(i_int) - t_interval(1);
pbs = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pbs = pbs(i_int);
pbs(1:safezone) = pbs(safezone);
pbd = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pbd = pbd(i_int);
pbd(1:safezone) = pbd(safezone);
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
pbm = pbm(i_int);
pbm(1:safezone) = pbm(safezone);
hr = dymget(dl, 'HR')/bpm2SI;
hr = hr(i_int);
hr(1:safezone) = hr(safezone);
co = dymget(dl, 'CO')*60000;
co = co(i_int);
co(1:safezone) = co(safezone);
lsv = dymget(dl, 'SV')/ml2SI;
tilt_a = dymget(dl, 'SystemicComponent.Tilt');
i_hut = find(tilt_a > 0, 1)
i_hutStop = find(tilt_a == max(tilt_a), 1)
% rsv = dymget(dl, 'heartComponent.pulmonaryValve.SV')/ml2SI;

psv = dymget(dl, 'P_sv')/mmHg2SI;
psv = psv(i_int);
ppv = dymget(dl, 'P_pv')/mmHg2SI;
ppv = ppv(i_int);

%% 
fig = figure(2);clf;
fs = 8;

s1 = subplot(4, 1, 1); cla; hold on;
set(gca, 'FontSize', fs);
title('60° HUT maneuver', 'FontSize', fs + 2)

pbp = fill([time;flipud(time)], [pbs;flipud(pbd)], color_lb,'EdgeColor',color_lb);
pbpm = plot(time, pbm, 'Color',color_b, 'LineWidth', 2);
% phr = plot(t, hr, 'r', 'LineWidth', 1);

% rescale the data
sc = 2/3;
sh = 45;
pas_d = T.PAs*sc + sh;
pad_d = T.PAd*sc + sh;
pbs_data = plot(T.PAs_t(i_past:end) + 10, pas_d(i_past:end), 'Color',color_s, 'LineWidth', 1);
pbd_data = plot(T.PAd_t(i_padt:end) + 10, pad_d(i_padt:end), 'Color',color_s, 'LineWidth', 1);

leg = legend([pbp(1), pbpm, pbs_data], 'model', 'mean (model)', 'data', 'Location', 'SouthEast', 'Orientation', 'Horizontal');
% leg = legend([pbp(1), pbpm, phr], 'PA [mmHg]', 'PA mean [mmHg]','HR [BPM]', 'Location', 'NorthEast');
leg.ItemTokenSize = [10, 2];
set(gca,'xtick',[]);
xlim(t_interval);
ylim([60, 130]);
yl = ylabel('PA (mmHg)');
yl_pos=get(yl,'Pos')
yl_x = yl_pos(1);

s1.Clipping = 'off';


s2 = subplot(4, 1, 2); hold on;
phr = plot(time, hr, 'Color',color_r , 'LineWidth', 1);
phr_data = plot(T.HR_t(i_hrt:end) + 10, T.HR(i_hrt:end), 'Color',color_s , 'LineWidth', 1);
xlim(t_interval);
ylim([60, 100]);
yl = ylabel('HR (bpm)');
yl_pos=get(yl,'Pos');
set(yl,'Pos',[yl_x yl_pos(2) yl_pos(3)]);
leg = legend('model', 'data', 'orientation', 'horizontal');
leg.ItemTokenSize = [10, 2];
set(gca,'xtick',[]);
set(gca,'ytick',[70, 80, 90]);

s3 = subplot(4, 1, 3); hold on;
set(gca, 'FontSize', fs);
plot(time, co, 'Color',color_m, 'LineWidth', 1)
co_base = co(find(co>1, 1)) + 0.4;
COf_data = plot(T.COf_t(i_coft:end) + 10, T.COf(i_coft:end)*co_base/100, 'Color',color_s , 'LineWidth', 1);
% plot(t(i_int), rsv(i_int), 'r')
set(gca,'xtick',[])
leg = legend('model', 'data', 'Location', 'NorthEast', 'orientation', 'horizontal')
leg.ItemTokenSize = [10, 2];
xlim(t_interval);
ylim([4, 8]);
yl = ylabel('CO (L/min)');
yl_pos=get(yl,'Pos');
set(yl,'Pos',[yl_x yl_pos(2) yl_pos(3)]);
set(gca,'xtick',[]);
s2.Clipping = 'off';

s4 = subplot(4, 1, 4); hold on;
set(gca, 'FontSize', fs);
% plot(t, tp, 'b')
% plot(t, pdv, 'r')
fsv = plot(time, psv, 'Color',color_b, 'LineWidth', 1);
fpv = plot(time, ppv, 'Color',color_r, 'LineWidth', 1);
xlim(t_interval);
yl = ylabel('P (mmHg)');
yl_pos=get(yl,'Pos')
set(yl,'Pos',[yl_x yl_pos(2) yl_pos(3)]);

xlabel('t (s)');
plot([time(i_hut), time(i_hut)], [50, 0], 'k--', 'LineWidth', 1)
plot([time(i_hutStop), time(i_hutStop)], [50, 0], 'k--', 'LineWidth', 1)
ylim([0, 12])
s4.Clipping = 'off';
leg = legend([fsv, fpv], 'Systemic veins', 'Pulmonary veins', 'Location', 'NorthEast', 'Orientation', 'horizontal');
leg.ItemTokenSize = [10, 2];


drawnow()
% move around
%%
vm = 0.8;
hm = 1;
tw = 7;
th = 11;
xax = 0.4;
hs = tw-1.2*hm;
vs = 0.25*th - vm - 0.1 + xax;
set(fig, 'Units', 'Centimeters', 'Position', [2, 4, tw, th])
set(s1, 'Units', 'Centimeters', 'Position', [hm, 0.75*th - xax, hs,vs + 0.3]);
set(s2, 'Units', 'Centimeters', 'Position', [hm, 0.5*th, hs,vs]);
set(s3, 'Units', 'Centimeters', 'Position', [hm, 0.25*th + xax, hs,vs]);
set(s4, 'Units', 'Centimeters', 'Position', [hm, vm, hs,vs]);
drawnow()
%%
s1.YLabel.FontSize = fs;

s3.YLabel.FontSize = fs;
s4.YLabel.FontSize = fs;
s2.YLabel.FontSize = fs;
s2.YAxis.FontSize = 8;
s2.YLabel.FontSize = fs;
s1.Legend.FontSize = fs;
s2.Legend.FontSize = fs;
s3.Legend.FontSize = fs;
s4.Legend.FontSize = fs;
%%

%% save
exportgraphics(gcf,'fig_R_Tilt.png','Resolution',300)
exportgraphics(gcf,'fig_R_Tilt.pdf', 'ContentType','vector')