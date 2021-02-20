%% plots the valsalva results
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')

datafile = '../../Results2/imp_tiltBase.mat';
dl = dymload(datafile)
%%
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
%%
time = dymget(dl, 'Time');
t_interval = [0, 60]; % interval in seconds
td = t_interval(2) - t_interval(1);
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int = time >= 2; % get rid of the initial zeros
safezone = 400;
i_hut = find(time > 20, 1)
i_hutStop = find(time > 21, 1)

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
% rsv = dymget(dl, 'heartComponent.pulmonaryValve.SV')/ml2SI;

psv = dymget(dl, 'P_sv')/mmHg2SI;
ppv = dymget(dl, 'P_pv')/mmHg2SI;

%% 
fig = figure(2);clf;
fs = 8;

s1 = subplot(2, 1, 1); cla; hold on;
set(gca, 'FontSize', fs);
title('60° HUT maneuver', 'FontSize', fs + 2)

pbp = fill([t;flipud(t)], [pbs;flipud(pbd)], [182/255, 226/255, 1],'EdgeColor',[182/255, 226/255, 1]);
pbpm = plot(t(i_int), pbm(i_int), 'b', 'LineWidth', 2);
phr = plot(t, hr, 'r', 'LineWidth', 1);
leg = legend([pbp(1), pbpm, phr], 'PA [mmHg]', 'PA mean [mmHg]','HR [BPM]', 'Location', 'NorthEast');
leg.ItemTokenSize = [10, 2];
set(gca,'xtick',[])
xlim([0,120])
ylim([60, 130])
s1.Clipping = 'off';

s2 = subplot(4, 1, 3); hold on;
set(gca, 'FontSize', fs);
plot(t(i_int), co(i_int), 'm', 'LineWidth', 1)
% plot(t(i_int), rsv(i_int), 'r')
set(gca,'xtick',[])
leg = legend('CO', 'Location', 'NorthEast')
leg.ItemTokenSize = [10, 2];
xlim([0,120])
ylabel('Flow L/min')
s2.Clipping = 'off';

s3 = subplot(4, 1, 4); hold on;
set(gca, 'FontSize', fs);
% plot(t, tp, 'b')
% plot(t, pdv, 'r')
fsv = plot(t, psv, 'b', 'LineWidth', 1)
fpv = plot(t, ppv, 'r', 'LineWidth', 1)
xlim([0,120])
ylabel('Pressure (mmHg)')
xlabel('t (s)')
plot([t(i_hut), t(i_hut)], [65, 0], 'k--', 'LineWidth', 1)
plot([t(i_hutStop), t(i_hutStop)], [65, 0], 'k--', 'LineWidth', 1)
ylim([0, 12])
s3.Clipping = 'off';
leg = legend([fsv, fpv], 'VC', 'PV', 'Location', 'NorthEast');
leg.ItemTokenSize = [10, 2];


drawnow()
% move around
vm = 0.8;
hm = 1;
tw = 7;
xax = 0.4;
set(fig, 'Units', 'Centimeters', 'Position', [2, 4, tw, th])
set(s1, 'Units', 'Centimeters', 'Position', [hm, 0.5*th, tw-1.2*hm,0.5*th - 0.5*vm]);
set(s2, 'Units', 'Centimeters', 'Position', [hm, 0.25*th + xax, tw-1.2*hm,0.25*th - 0.1 - xax]);
set(s3, 'Units', 'Centimeters', 'Position', [hm, vm, tw-1.2*hm,0.25*th - vm - 0.1 + xax]);

%% save
exportgraphics(gcf,'fig_R_Tilt.png','Resolution',300)
exportgraphics(gcf,'fig_R_Tilt.pdf', 'ContentType','vector')