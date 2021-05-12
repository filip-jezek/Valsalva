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
% read and plot valsalva data from Kosinksi et al. 2018
data1 = xlsread('..\..\data\Valsalva\norm004_BP.xlsx', 'A3:C98711');
data2 = xlsread('..\..\data\Valsalva\norm004_HR.xlsx', 'A5:M580');
t3 = xlsread('..\..\data\Valsalva\norm004_Pexp.xlsx', 'B2:JN2');
%%
data_interval1 = 254;
data_interval2 = data_interval1 + 60;
% tps = xlsread('..\..\data\Valsalva\norm004_Pexp.xlsx', 'B6:JN15');
data1_i = data1(:, 1) >= data_interval1 & data1(:, 1) < data_interval2;
data2_i = data2(:, 1) >= data_interval1 & data2(:, 1) < data_interval2;
t1_vm = data1(data1_i, 1) - data_interval1;
% finap = data1(:, 2);
reBap = data1(data1_i, 3);
t2_vm = data2(data2_i, 1)- data_interval1;
hr_vm = data2(data2_i, 8);
% lvet = data2(data2_i, 10);
sv_vm = data2(data2_i, 11);
% tpr = data2(data2_i, 12);
% tp = mean(tps, 1);

% figure(1);clf;hold on;
% plot(t1_vm, reBap);
% plot(t2_vm, hr_vm);
% % plot(t2_vm, lvet);
% plot(t2_vm, sv_vm);
% % plot(t2_vm, tpr);
% [vmkos_upenv,vmkos_lowend] = envelope(reBap, 100, 'peak');
% % plot(t1_vm, vmkos_upenv, t1_vm, vmkos_lowend);
% fill([t1_vm; flipud(t1_vm)], [vmkos_upenv; flipud(vmkos_lowend(201:end));repmat(vmkos_lowend(200), [200, 1])], [0.8, 0.8, 0.8],'EdgeColor',[182/255, 226/255, 1])

%%
fig = figure(2);clf;
fs = 8;
%%
s1 = subplot(2, 1, 1); hold on;
set(gca, 'FontSize', fs);
title('Valsalva maneuver', 'FontSize', fs + 2)

% draw data for comparison
[vmkos_upenv,vmkos_lowend] = envelope(reBap, 100, 'peak');
% p_vm = plot(t1_vm, vmkos_upenv, t1_vm, vmkos_lowend, 'Color', [0.8, 0.8, 0.8]);
p_vm = fill([t1_vm; flipud(t1_vm)], [vmkos_upenv; flipud(vmkos_lowend(201:end));repmat(vmkos_lowend(200), [200, 1])], [0.8, 0.8, 0.8],'EdgeColor',[0.8, 0.8, 0.8]);
vm_hr = plot(t2_vm, hr_vm, 'k', 'Linewidth', 1);
% draw the simulation output
pbp = plot(t, pb, 'b');
pbpm = plot(t(i_int), pbm(i_int), 'b', 'LineWidth', 2);
phr = plot(t, hr, 'r', 'LineWidth', 1);
plot(t1_vm, vmkos_upenv, t1_vm, vmkos_lowend);
leg = legend([p_vm, pbp, vm_hr, phr], 'PA [mmHg] (data)', 'PA [mmHg] (model)', 'HR [BPM] (data)', 'HR [BPM] (model)', 'Location', 'NorthWest');
leg.ItemTokenSize = [10, 10];
ylim([45, 160])
set(gca,'xtick',[])
s1.Clipping = 'off';
%%
s2 = subplot(4, 1, 3); cla;hold on;
set(gca, 'FontSize', fs);
plot(t2_vm, sv_vm, 'Color', [0.8, 0.8, 0.8], 'Linewidth', 2);
plot(t(i_int), lsv(i_int), 'b', 'LineWidth', 1)
plot(t(i_int), rsv(i_int), 'r', 'LineWidth', 1)
set(gca,'xtick',[], 'ytick', [50, 100, 150, 200])
leg = legend('LV SV [data]', 'LV SV (model)', 'RV SV (model)', 'Location', 'NorthWest');
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