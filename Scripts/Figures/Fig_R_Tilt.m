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
time = dymget(dl, 'Time');
t_interval = [0, 300]; % interval in seconds
td = t_interval(2) - t_interval(1);
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int = time >= 2; % get rid of the initial zeros
safezone = 400;


t = time - t_interval(1);
pbs = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pbs(1:safezone) = pbs(safezone);
pbd = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pbd(1:safezone) = pbd(safezone);
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
pbm(1:safezone) = pbm(safezone);
hr = dymget(dl, 'HR')/bpm2SI;
co = dymget(dl, 'CO')*60000;
lsv = dymget(dl, 'SV')/ml2SI;
tilt_a = dymget(dl, 'SystemicComponent.Tilt');
i_hut = find(tilt_a > 0, 1)
i_hutStop = find(tilt_a == max(tilt_a), 1)
% rsv = dymget(dl, 'heartComponent.pulmonaryValve.SV')/ml2SI;

psv = dymget(dl, 'P_sv')/mmHg2SI;
ppv = dymget(dl, 'P_pv')/mmHg2SI;

%% 
fig = figure(2);clf;
fs = 8;

s1 = subplot(4, 1, 1); cla; hold on;
set(gca, 'FontSize', fs);
title('60° HUT maneuver', 'FontSize', fs + 2)

pbp = fill([t;flipud(t)], [pbs;flipud(pbd)], color_lb,'EdgeColor',color_lb);
pbpm = plot(t(i_int), pbm(i_int), 'Color',color_b, 'LineWidth', 2);
% phr = plot(t, hr, 'r', 'LineWidth', 1);
leg = legend([pbp(1), pbpm], 'PA', 'PA mean', 'Location', 'NorthEast', 'Orientation', 'Horizontal');
% leg = legend([pbp(1), pbpm, phr], 'PA [mmHg]', 'PA mean [mmHg]','HR [BPM]', 'Location', 'NorthEast');
leg.ItemTokenSize = [10, 2];
set(gca,'xtick',[]);
xlim(t_interval);
ylim([80, 130]);
yl = ylabel('PA (mmHg)');
yl_pos=get(yl,'Pos')
yl_x = yl_pos(1);

s1.Clipping = 'off';


s2 = subplot(4, 1, 2); hold on;
phr = plot(t, hr, 'Color',color_r , 'LineWidth', 1);
xlim(t_interval);
ylim([60, 100]);
yl = ylabel('HR (bpm)');
yl_pos=get(yl,'Pos');
set(yl,'Pos',[yl_x yl_pos(2) yl_pos(3)]);
leg = legend('HR');
leg.ItemTokenSize = [10, 2];
set(gca,'xtick',[]);
set(gca,'ytick',[70, 80, 90]);

s3 = subplot(4, 1, 3); hold on;
set(gca, 'FontSize', fs);
plot(t(i_int), co(i_int), 'Color',color_m, 'LineWidth', 1)
% plot(t(i_int), rsv(i_int), 'r')
set(gca,'xtick',[])
leg = legend('CO', 'Location', 'NorthEast')
leg.ItemTokenSize = [10, 2];
xlim(t_interval);
yl = ylabel('CO (L/min)');
yl_pos=get(yl,'Pos');
set(yl,'Pos',[yl_x yl_pos(2) yl_pos(3)]);
set(gca,'xtick',[]);
s2.Clipping = 'off';

s4 = subplot(4, 1, 4); hold on;
set(gca, 'FontSize', fs);
% plot(t, tp, 'b')
% plot(t, pdv, 'r')
fsv = plot(t, psv, 'Color',color_b, 'LineWidth', 1);
fpv = plot(t, ppv, 'Color',color_r, 'LineWidth', 1);
xlim(t_interval);
yl = ylabel('P (mmHg)');
yl_pos=get(yl,'Pos')
set(yl,'Pos',[yl_x yl_pos(2) yl_pos(3)]);

xlabel('t (s)');
plot([t(i_hut), t(i_hut)], [50, 0], 'k--', 'LineWidth', 1)
plot([t(i_hutStop), t(i_hutStop)], [50, 0], 'k--', 'LineWidth', 1)
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