%% plots all the figures for figure summary
%% color defiinitons
color_b = [28, 108, 200]/255;
color_r = [238, 46, 47]/255;
color_g = [0, 140, 72]/255;
color_m = [226, 113, 199]/255;
color_lb = [182, 226, 255]/255;
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
%% plots the baseline results
datafile = '../../Results/CardiovascularSystem.mat'
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
dl = dymload(datafile)

time_bl = dymget(dl, 'Time');
t_interval_bl = [28.3, 29.3]; % interval in seconds
t_interval_bl2 = [27.95, 29.95]; % interval in seconds
td = t_interval_bl(2) - t_interval_bl(1)
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int_bl = (time_bl >= t_interval_bl(1) & time_bl <= t_interval_bl(2));

t_bl = time_bl - t_interval_bl(1);
A_p1 = dymget(dl, 'SystemicComponent.ascending_aorta_A.p1')/mmHg2SI;
A_p2 = dymget(dl, 'SystemicComponent.internal_carotid_R8_B.p_C')/mmHg2SI;
A_p3 = dymget(dl, 'brachial_pressure')/mmHg2SI;

B1_vlv = dymget(dl, 'V_LV')/ml2SI;
B1_vla = dymget(dl, 'V_la')/ml2SI;

B1_pv = dymget(dl, 'heartComponent.pv.pressure')/mmHg2SI;
B1_psa = A_p1;
B1_plv = dymget(dl, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
B1_pla = dymget(dl, 'heartComponent.mitralValve.q_in.pressure')/mmHg2SI;

B1_vrv = dymget(dl, 'heartComponent.ventricles.V_RV')/ml2SI;
B1_prv = dymget(dl, 'heartComponent.ventricles.P_RV')/mmHg2SI;
B1_vra = dymget(dl, 'heartComponent.ra.volume')/ml2SI;


psa = dymget(dl, 'heartComponent.sa.pressure')/mmHg2SI;
plv = dymget(dl, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
pla = dymget(dl, 'heartComponent.mitralValve.q_in.pressure')/mmHg2SI;

vlv = dymget(dl, 'V_LV')/ml2SI;
vla = dymget(dl, 'V_la')/ml2SI;

ppa = dymget(dl, 'heartComponent.pa.pressure')/mmHg2SI;
prv = dymget(dl, 'heartComponent.tricuspidValve.q_out.pressure')/mmHg2SI;
pra = dymget(dl, 'heartComponent.tricuspidValve.q_in.pressure')/mmHg2SI;

vrv = dymget(dl, 'heartComponent.ventricles.V_RV')/ml2SI;
vra = dymget(dl, 'heartComponent.ra.volume')/ml2SI;

pa_tibial = dymget(dl, 'SystemicComponent.posterior_tibial_T4_L214.p_in')/mmHg2SI;
pv_tibial = dymget(dl, 'SystemicComponent.posterior_tibial_T4_L214.p_out')/mmHg2SI;
qmv = dymget(dl, 'q_mv')/ml2SI/1000*60;
%% SIZE
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 5.8;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])
hold on;
%% A1
% title('A: Left ventricular pressures and volumes', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% title('A: Left ventricular pressures and volumes');
plot(t_bl(i_int_bl), A_p1(i_int_bl), 'LineWidth', 1, 'Color', color_r);
plot(t_bl(i_int_bl), A_p2(i_int_bl), 'LineWidth', 1, 'Color', color_b);
plot(t_bl(i_int_bl), A_p3(i_int_bl), 'LineWidth', 1, 'Color', color_g);
% plot(t(i_int), pb(i_int));
% set(gca,'xtick',[])
leg = legend('Asc Aor', 'Int carotid', 'Brach art', 'Location', 'NorthEast')
leg.ItemTokenSize = [10, 2];
% legend('boxoff')
ylim([80, 130]);
ylabel('Pressure (mmHg)')
xlabel('t (s)')
% s_a1.Clipping = 'off';
xlim([0 td])

exportgraphics(gcf,'fig_R_sum_A1.png','Resolution',300)
exportgraphics(gcf,'fig_R_SUM_A1.pdf', 'ContentType','vector')

%% A2
figure(1);clf; hold on;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 6.2;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])

plot(t_bl(i_int_bl), B1_pv(i_int_bl), 'LineWidth', 1, 'Color', color_m);
plot(t_bl(i_int_bl), B1_pla(i_int_bl), 'LineWidth', 1, 'Color', color_g);
plot(t_bl(i_int_bl), B1_plv(i_int_bl), 'LineWidth', 1, 'Color', color_r);
plot(t_bl(i_int_bl), B1_psa(i_int_bl), 'LineWidth', 1, 'Color', color_b);
ylim([0, 120]);
ylabel('Pressure (mmHg)')
xlabel('t (s)')
yyaxis right;cla;
% set(gca,'YColor','k');
plot(t_bl(i_int_bl), B1_vla(i_int_bl), ':', 'LineWidth', 1, 'Color', color_g);
plot(t_bl(i_int_bl), B1_vlv(i_int_bl), ':', 'LineWidth', 1, 'Color', color_r);
g = gca; g.YAxis(2).Color = [0,0,0];
ylim([0, 160]);
yticks([0, 50, 100])
ylabel('Volume (mL)');

leg = legend('P PV', 'P LA', 'P LV','P Asc Aor', 'V LA*', 'V LV*', 'Location', 'East');
leg.ItemTokenSize = [10, 2];
% legend('boxoff')
xlim([0 td])
exportgraphics(gcf,'fig_R_SUM_A2.pdf', 'ContentType','vector')

%% A3
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 5.8;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])
hold on;

% volumes
plot(B1_vlv(i_int_bl), B1_plv(i_int_bl), 'LineWidth', 1, 'Color', color_b);
plot(B1_vrv(i_int_bl), B1_prv(i_int_bl), 'LineWidth', 1, 'Color', color_r);

% plot(vla(i_int_bl), pla(i_int_bl), 'LineWidth', 1, 'Color', color_g);
% plot(vra(i_int_bl), pra(i_int_bl), 'LineWidth', 1, 'Color', color_m);
ylim([0, 120]);
xlim([40 160])
ylabel('Pressure (mmHg)');
xlabel('Volume (mL)')
leg = legend('LV', 'RV', 'Location', 'best');
leg.ItemTokenSize = [10, 2];

exportgraphics(gcf,'fig_R_SUM_A3.pdf', 'ContentType','vector')
%% B - tilt data read
datafile = '../../Results/CVS_tiltable.mat';
% system('"c:/Program Files/Dymola 2021x/bin/dsres2sdf.exe" ../../Results/CVS_tiltable.mat ../../Results/CVS_tiltable.sdf');
% dl_arst = '../../Results/CVS_tiltable.sdf'

% time_arst = h5read(dl_arst, '/Time');
% %%
dl = dymload(datafile)
time = dymget(dl, 'Time');
t_interval_hut = [0, 600]; % interval in seconds
td = t_interval_bl(2) - t_interval_bl(1);
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int_hut = time >= 2; % get rid of the initial zeros
safezone = 400;

t_hut = time - t_interval_hut(1);
pbs = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pbs(1:safezone) = pbs(safezone);
pbd = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pbd(1:safezone) = pbd(safezone);
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
pbm(1:safezone) = pbm(safezone);
pa_tibial_hut = dymget(dl, 'SystemicComponent.posterior_tibial_T4_L214.p_in')/mmHg2SI;
pv_tibial_hut = dymget(dl, 'SystemicComponent.posterior_tibial_T4_L214.p_out')/mmHg2SI;

hr = dymget(dl, 'HR')/bpm2SI;
co = dymget(dl, 'CO')*60000;
lsv = dymget(dl, 'SV')/ml2SI;
% rsv = dymget(dl, 'heartComponent.pulmonaryValve.SV')/ml2SI;

%% B1 - tilt data plot
figure(1);clf; hold on;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 6.2;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])

psv = dymget(dl, 'P_sv')/mmHg2SI;
ppv = dymget(dl, 'P_pv')/mmHg2SI;

pbp = fill([t_hut;flipud(t_hut)]/60, [pbs;flipud(pbd)], color_lb,'EdgeColor',[182/255, 226/255, 1]);
pbpm = plot(t_hut(i_int_hut)/60, pbm(i_int_hut), 'LineWidth', 2, 'Color', color_b);
ylabel('Pressure (mmHg)');
ylim([60, 130]);

yyaxis right;
% set(gca,{'ycolor'},{'r';'r'})
phr = plot(t_hut/60, hr, 'r', 'LineWidth', 1, 'Color', color_r);
% pco = plot(t(i_int)/60, co(i_int), 'LineWidth', 1, 'Color', color_m);
ylim([60, 130])
g = gca; g.YAxis(2).Color = [0,0,0];
ylabel('HR (bpm)');
xlabel('t (min)');
leg = legend('PA', 'PA mean', 'HR*', 'Location', 'SouthEast', 'orientation', 'horizontal');
leg.ItemTokenSize = [10, 2];

exportgraphics(gcf,'fig_R_SUM_B1.pdf', 'ContentType','vector')
%% B2 - tilt peripherals
figure(1);clf; hold on;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 5.8;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])

b2_d = 2
i_bl = find(time_bl >= 27.8, 1)
i_bl_stop = find(time_bl >= 27.8 + b2_d, 1)

i_hut = find(time > 120, 1);
i_hutStop = find(time > 120 + b2_d, 1);

int_bl = [i_bl : i_bl_stop];
int_hut = [i_hut:i_hutStop];
plot(time_bl(int_bl) - time_bl(i_bl), pa_tibial(int_bl), ':', 'Color', color_r);
plot(time_bl(int_bl) - time_bl(i_bl), pv_tibial(int_bl), ':', 'Color', color_b);
plot(t_hut(int_hut) - time(i_hut), pa_tibial_hut(int_hut), 'Color', color_r);
plot(t_hut(int_hut)- time(i_hut), pv_tibial_hut(int_hut), 'Color', color_b);
ylim([0, 220])
leg = legend('Art. (sup)', 'Ven. (sup)', 'Art. (HUT)', 'Ven. (HUT)', 'Location', 'SouthEast', 'Orientation', 'vertical');
leg.ItemTokenSize = [10, 2];
ylabel('Pressure (mmHg)');
xlabel('t (s)');
exportgraphics(gcf,'fig_R_SUM_B2.pdf', 'ContentType','vector')

%% C - valsalva
datafile = '../../Results/CVS_valsalva.mat';
dl = dymload(datafile);

time = dymget(dl, 'Time');
t_interval = [0, 60]; % interval in seconds
td = t_interval(2) - t_interval(1);
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int = time >= 2; % get rid of the initial zeros

t = time - t_interval(1);
pb = dymget(dl, 'brachial_pressure')/mmHg2SI;
pbs = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pbd = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
hr = dymget(dl, 'HR')/bpm2SI;
% lsv = dymget(dl, 'SV')/ml2SI;
% rsv = dymget(dl, 'heartComponent.pulmonaryValve.CO')./dymget(dl, 'HR')/ml2SI;

% plot C
figure(1);clf; hold on;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 5.8;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])
x = 400;
p_vm = fill([time(x:end); flipud(time(x:end))], [pbs(x:end); flipud(pbd(x:end))], color_lb,'EdgeColor',color_lb);
% draw the simulation output
% pbp = plot(t_bl, pb, 'b');
pbpm = plot(t(i_int), pbm(i_int), 'LineWidth', 1, 'Color', color_b);
phr = plot(t, hr, 'LineWidth', 1, 'Color', color_r);

leg = legend('PA (mmHg)', 'PA mean (mmHg)', 'HR (bpm)', 'Location', 'NorthWest');
leg.ItemTokenSize = [10, 10];
ylim([50, 170])
yticks([60, 100, 140])
ylabel(' ');
xlabel('t (s)')
exportgraphics(gcf,'fig_R_SUM_C.pdf', 'ContentType','vector');

%% D: hemorrhage
datafile = '../../Results/Hemorrhage.mat';
dl = dymload(datafile);

time = dymget(dl, 'Time');
t_interval = [0, time(end)]; % interval in seconds
td = t_interval(2) - t_interval(1);
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int = time >= 2; % get rid of the initial zeros

t = time - t_interval(1);
pb = dymget(dl, 'brachial_pressure')/mmHg2SI;
pbs = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pbd = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
co = dymget(dl, 'CO')*60000;
hr = dymget(dl, 'HR')/bpm2SI;
sprintf('Hemorrhage at %d mins: HR %d, BPm %d, CO %d lpm \n', time(end) /60,hr(end), pbm(end), co(end))

figure(1);clf; hold on;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 5.8;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])
x = 20;
p_m = fill([time(x:end); flipud(time(x:end))]/60, [pbs(x:end); flipud(pbd(x:end))], color_lb,'EdgeColor',color_lb);
% draw the simulation output
% pbp = plot(t_bl, pb, 'b');
pbpm = plot(t(x:end)/60, pbm(x:end), 'LineWidth', 1, 'Color', color_b);
phr = plot(t(x:end)/60, hr(x:end), 'LineWidth', 1, 'Color', color_r);

leg = legend('PA (mmHg)', 'PA mean (mmHg)', 'HR (bpm)', 'Location', 'northeast');
leg.ItemTokenSize = [10, 10];
ylim([60, 140])
yticks([60, 80, 100, 120])
ylabel(' ');
xlabel('t (min)')
exportgraphics(gcf,'fig_R_SUM_D.pdf', 'ContentType','vector');

%% E: Exercise
datafile = '../../Results/ExStep_base.mat';
dl  = dymload(datafile)
time = dymget(dl, 'Time');

pb = dymget(dl, 'brachial_pressure')/mmHg2SI;
pbs = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pbd = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pbm = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
hr = dymget(dl, 'HR')/bpm2SI;
SV = dymget(dl, 'SV')/ml2SI;
CO = dymget(dl, 'CO')*60000;
el = dymget(dl, 'Exercise.y_');

ts = 35
te = 36.6
td = te - ts;
diffs = [1; find( diff(el) > 0)];
ef = [];
edv = [];
pow_lvs = [];
pow_rvs = [];
co_s = [];
q_ex_s = [];
sv_s = [];

pbs_s = [];
pbd_s = [];
pbm_s = [];

for i = 1:size(diffs)
    i_s = find(time > time(diffs(i)) + ts, 1);
    i_e = find(time > time(diffs(i)) + te, 1);
    i_int = [i_s:i_e];
    co_s = [co_s, mean(CO(i_int))];
    sv_s = [sv_s, mean(SV(i_int))];
    pbs_s = [pbs_s, mean(pbs(i_int))];
    pbd_s = [pbd_s, mean(pbd(i_int))];
    pbm_s = [pbm_s, mean(pbm(i_int))];
end

figure(1);clf; hold on;
set(gcf, 'DefaultAxesFontSize', 6, 'defaultLineLineWidth',1.0);
tw = 6.2;
th = 4.5;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th])

% plot([0:10:100], edv, 'b*-')
% plot([0:10:100], pow_rvs, 'r*-')
% plot([0:10:100], pow_lvs, 'b*-')
% plot([0:10:100], ef*10, 'b*-')
el_s = [0:10:100]';
fill([el_s; flipud(el_s)], [pbs_s'; flipud(pbd_s')], color_lb,'EdgeColor',color_lb);
plot(el_s, pbm_s, 'Color', color_b);
plot(el_s, sv_s, 'Color', color_g);
ylabel('Pressure (mmHg)')
ylim([50, 250])
yyaxis right;
plot([0:10:100], co_s, '--', 'Color', color_m, 'LineWidth', 1.5);
ylim([5, 25])
ylabel('CO (L/min)')
xlabel('Exercise level (% of max)')
g = gca; g.YAxis(2).Color = color_m;
leg = legend('PA', 'PA mean', 'CO*', 'Location', 'NorthWest');
leg.ItemTokenSize = [10, 10];
exportgraphics(gcf,'fig_R_SUM_E.pdf', 'ContentType','vector');

% plot([0:10:100], q_ex_s, 'r*-')
% title('Ejection fraction on exercise level')