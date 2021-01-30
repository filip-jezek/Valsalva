%% plots the valsalva results
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')

datafile = '../../Results2/imp_stepEx_normal.sdf';
hi = h5info(datafile)
%%
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
%%

dl = datafile
time = h5read(dl, '/Time');
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
%%
t = time;
pb = h5read(dl, '/brachial_pressure')/mmHg2SI;

pbs = h5read(dl, '/brachial_pressure_systolic')/mmHg2SI;
pbs(1:safezone) = pbs(safezone);
pbd = h5read(dl, '/brachial_pressure_diastolic')/mmHg2SI;
pbd(1:safezone) = pbd(safezone);
pbm = h5read(dl, '/brachial_pressure_mean')/mmHg2SI;
pbm(1:safezone) = pbm(safezone);
psa = h5read(dl, '/heartComponent/sa/pressure')/mmHg2SI;

hr = h5read(dl, '/HR')/bpm2SI;
lsv = h5read(dl, '/SV')/ml2SI;
% rsv = h5read(dl, '/heartComponent/pulmonaryValve/SV')/ml2SI;
plv = h5read(dl, '/P_LV')/mmHg2SI;
prv = h5read(dl, '/heartComponent/ventricles/P_RV')/mmHg2SI;
vlv = h5read(dl, '/heartComponent/ventricles/V_LV')/ml2SI;
vrv = h5read(dl, '/heartComponent/ventricles/V_RV')/ml2SI;
vla = h5read(dl, '/V_la')/ml2SI;

pla = h5read(dl, '/P_LA')/mmHg2SI;

psv = h5read(dl, '/P_sv')/mmHg2SI;
ppv = h5read(dl, '/P_pv')/mmHg2SI;
ppa = h5read(dl, '/heartComponent/pa/pressure')/mmHg2SI;

pra = h5read(dl, '/heartComponent/tricuspidValve/q_in/pressure')/mmHg2SI;
vra = h5read(dl, '/heartComponent/ra/volume')/ml2SI;


tp = h5read(dl, '/thoracic_pressure')/mmHg2SI;
pdv = h5read(dl, '/SystemicComponent/femoral_vein_R34/port_a/pressure')/mmHg2SI;
psv = h5read(dl, '/P_sv')/mmHg2SI - tp;
ppv = h5read(dl, '/P_pv')/mmHg2SI - tp;
el = h5read(dl, '/Exercise/y_')*100;
qmv = h5read(dl, '/q_mv')/ml2SI/1000*60;
%% plot pv loops
ts = 35
te = 36
diffs = [1; find( diff(el) > 0)]
figure(3); clf;hold on;
ef = [];
edv = [];
for i = 1:size(diffs)
%     i = 9

    i_s = find(time > time(diffs(i)) + ts, 1);
    i_e = find(time > time(diffs(i)) + te, 1);
    t = time - time(i_s);
i_int = [i_s:i_e];
%     plot(vlv(i_s:i_e), plv(i_s:i_e))
%     plot(vrv(i_s:i_e), prv(i_s:i_e))
    plot(vla(i_s:i_e), pla(i_s:i_e))
    sv =  max(vlv(i_s:i_e)) -  min(vlv(i_s:i_e));
    ef = [ef, sv / max(vlv(i_s:i_e))];
    edv = [edv, min(vlv(i_s:i_e))];
    
end

plot([0:10:100], edv, 'b*-')
% plot([0:10:100], ef*100, 'b*-')
% title('Ejection fraction on exercise level')
%%

fig = figure(2);clf;
fs = 8;

s1 = subplot(2, 1, 1); hold on;
set(gca, 'FontSize', fs);
title('Valsalva maneuver', 'FontSize', fs + 2)
pbp = plot(t, pb, 'b');
pbpm = plot(t(i_int), pbm(i_int), 'b', 'LineWidth', 2);
phr = plot(t, hr, 'r');
legend([pbp, phr], 'P BA [mmHg]', 'HR [BPM]', 'Location', 'NorthWest');
set(gca,'xtick',[])
s1.Clipping = 'off';

s2 = subplot(4, 1, 3); hold on;
set(gca, 'FontSize', fs);
plot(t(i_int), lsv(i_int), 'b')
plot(t(i_int), rsv(i_int), 'r')
set(gca,'xtick',[], 'ytick', [50, 100, 150])
legend('LV SV', 'RV SV', 'Location', 'NorthWest')
ylim([0, 200])
ylabel('ml')
s2.Clipping = 'off';

s3 = subplot(4, 1, 4); hold on;
set(gca, 'FontSize', fs);
plot(t, tp, 'b')
plot(t, pdv, 'r')
plot(t, psv, 'g')
plot(t, ppv, 'm')
ylabel('mmHg')
xlabel('time')
leg = legend('P Thor.', 'P DV', 'TP SV', 'TP PV', 'Location', 'NorthWest');
leg.ItemTokenSize = [10, 2];

drawnow()
%% move around
vm = 0.8;
hm = 1;
tw = 7;
th = 10;
set(fig, 'Units', 'Centimeters', 'Position', [2, 4, tw, th])
set(s1, 'Units', 'Centimeters', 'Position', [hm, 0.5*th, tw-1.2*hm,0.5*th - 0.5*vm]);
set(s2, 'Units', 'Centimeters', 'Position', [hm, 0.25*th, tw-1.2*hm,0.25*th - 0.1]);
set(s3, 'Units', 'Centimeters', 'Position', [hm, vm, tw-1.2*hm,0.25*th - vm - 0.1]);

%% save
exportgraphics(gcf,'fig_R_VM.png','Resolution',300)


%% wiggers
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 10);

s_a1 = subplot(4, 2, 1);
hold on;
title('A: LV pressures and volumes');
plot(t(i_int), plv(i_int), 'b');
plot(t(i_int), psa(i_int), 'r');
plot(t(i_int), pla(i_int), 'm');
plot(t(i_int), pb(i_int));
set(gca,'xtick',[])
legend('P LV', 'P Asc Aor', 'P brach art', 'P LA', 'Location', 'SouthEast')
% legend('boxoff')
ylim([0, 220]);
ylabel('Pressure [mmHg]')
s_a1.Clipping = 'off';
% xlim([0 td])
    
% volumes
s_a2 = subplot(4, 2, 3);hold on;
plot(t(i_int), vlv(i_int), 'b');
plot(t(i_int), vla(i_int), 'r');
set(gca,'xtickMode', 'auto')
ylim([0, 250]);
% xlim([0 td])
xlabel('Time [s]');
ylabel('Volume [ml]')
s_a1.Clipping = 'off';
% pos = get(s1, 'Position');
% set(s1, 'Position', [0.05, 0.5, 0.9,0.4])
% set(s2, 'Position', [0.05, 0.05, 0.9,0.4])

% B
s_b1 = subplot(4, 2, 2);cla;
hold on;
title('B: RV pressures and volumes');
plot(t(i_int), prv(i_int), 'b');
plot(t(i_int), ppa(i_int), 'r');
plot(t(i_int), pra(i_int), 'm');
% set(gca,'xtick',[], 'ytick', [])
set(gca,'xtick',[])
legend('P RV', 'P Pulm Ar', 'P RA')
% ylim([0, 120])
% xlim([0 td])
% volumes
s_b2 = subplot(4, 2, 4);hold on;
plot(t(i_int), vrv(i_int), 'b');
plot(t(i_int), vra(i_int), 'r');
set(gca,'xtickMode', 'auto', 'ytick', [])
ylim([0, 250]);
% xlim([0 td])
legend('V LV', 'V LA')
xlabel('Time [s]');
% pos = get(s1, 'Position');
% set(s1, 'Position', [0.05, 0.5, 0.9,0.4])
% set(s2, 'Position', [0.05, 0.05, 0.9,0.4])

% C
s_c = subplot(2, 2, 3);cla;hold on;
title('C: Mitral flow profile')
plot(t(i_int), qmv(i_int));
xlabel('Time [s]');
ylabel('Flow [L/min]');

% D
s_d = subplot(2, 2, 4);cla;hold on;
title('D: Ventricular PV loop')
% plot(vra(i_int), pra(i_int));
% plot(vla(i_int), pla(i_int));
plot(vrv(i_int), prv(i_int), 'r');
plot(vlv(i_int), plv(i_int), 'b');
xlabel('Time [s]');
ylabel('Pressure [mmHg]');

%% align the axis
drawnow()
pos = get(s_a2, 'Position');

hm = 0.1;
vm = 0.05;
set(s_a1, 'Position', [hm, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
set(s_a2, 'Position', [hm, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])

set(s_b1, 'Position', [0.5+hm/2, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
set(s_b2, 'Position', [0.5+hm/2, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])

set(s_c, 'Position', [hm, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])
set(s_d, 'Position', [0.5 + hm/2, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])
%%
tw = 17;
set(s_a1, 'Units', 'Centimeters', 'Position', tw*[hm, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
set(s_a2, 'Units', 'Centimeters', 'Position', tw*[hm, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])

set(s_b1, 'Units', 'Centimeters', 'Position', tw*[0.5+hm/2, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
set(s_b2, 'Units', 'Centimeters', 'Position', tw*[0.5+hm/2, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])

set(s_c, 'Units', 'Centimeters', 'Position', tw*[hm, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])
set(s_d, 'Units', 'Centimeters', 'Position', tw*[0.5 + hm/2, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])

%% Draw boxes
annotation(fig1,'rectangle', [0 0.5 0.5 0.5]);
annotation(fig1,'rectangle', [0.5 0.5 0.5 0.5]);


%% save figure
tw = 17;
th = 17;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th], 'PaperUnits', 'Centimeters', 'PaperSize', [tw, th])
print(gcf,'fig1_200.png', '-dpng', '-r200')
% saveas(gcf,'fig1.svg')