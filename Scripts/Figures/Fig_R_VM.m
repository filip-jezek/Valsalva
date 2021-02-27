%% plots the valsalva results
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')

datafile = '../../Results2/CVS_valsalva.mat';
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

t = time - t_interval(1);
pb = dymget(dl, 'brachial_pressure')/mmHg2SI;
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
hr = dymget(dl, 'HR')/bpm2SI;
lsv = dymget(dl, 'SV')/ml2SI;
rsv = dymget(dl, 'heartComponent.pulmonaryValve.CO')./dymget(dl, 'HR')/ml2SI;

tp = dymget(dl, 'thoracic_pressure')/mmHg2SI;
pdv = dymget(dl, 'SystemicComponent.femoral_vein_R34.port_a.pressure')/mmHg2SI;
psv = dymget(dl, 'P_sv')/mmHg2SI - tp;
ppv = dymget(dl, 'P_pv')/mmHg2SI - tp;

%% 
fig = figure(2);clf;
fs = 8;
%%
s1 = subplot(2, 1, 1); hold on;
set(gca, 'FontSize', fs);
title('Valsalva maneuver', 'FontSize', fs + 2)
pbp = plot(t, pb, 'b');
pbpm = plot(t(i_int), pbm(i_int), 'b', 'LineWidth', 2);
phr = plot(t, hr, 'r', 'LineWidth', 1);
leg = legend([pbp, phr], 'PA [mmHg]', 'HR [BPM]', 'Location', 'NorthWest');
leg.ItemTokenSize = [10, 10];
ylim([45, 160])
set(gca,'xtick',[])
s1.Clipping = 'off';

s2 = subplot(4, 1, 3); hold on;
set(gca, 'FontSize', fs);
plot(t(i_int), lsv(i_int), 'b', 'LineWidth', 1)
plot(t(i_int), rsv(i_int), 'r', 'LineWidth', 1)
set(gca,'xtick',[], 'ytick', [50, 100, 150, 200])
leg = legend('LV SV', 'RV SV', 'Location', 'NorthWest');
leg.ItemTokenSize = [10, 10];
ylim([0, 250])
ylabel('Volume (ml)')
s2.Clipping = 'off';

s3 = subplot(4, 1, 4); hold on;
set(gca, 'FontSize', fs);
plot(t, tp, 'b', 'LineWidth', 1)
plot(t, pdv, 'r', 'LineWidth', 1)
plot(t, psv, 'g', 'LineWidth', 1)
plot(t, ppv, 'm', 'LineWidth', 1)
ylabel('Pressure (mmHg)')
xlabel('t (s)')
leg = legend('P Thor.', 'P DV', 'TP SV', 'TP PV', 'Location', 'NorthWest');
leg.ItemTokenSize = [10, 2];

drawnow()
%% move around
vm = 0.8;
hm = 1;
tw = 7;
th = 11;
xax = 0.4;
set(fig, 'Units', 'Centimeters', 'Position', [2, 4, tw, th])
set(s1, 'Units', 'Centimeters', 'Position', [hm, 0.5*th, tw-1.2*hm,0.5*th - 0.5*vm]);
set(s2, 'Units', 'Centimeters', 'Position', [hm, 0.25*th + xax, tw-1.2*hm,0.25*th - 0.1 - xax]);
set(s3, 'Units', 'Centimeters', 'Position', [hm, vm, tw-1.2*hm,0.25*th - vm - 0.1 + xax]);

%% save
exportgraphics(gcf,'fig_R_VM.png','Resolution',300)
exportgraphics(gcf,'fig_R_VM.pdf', 'ContentType','vector')