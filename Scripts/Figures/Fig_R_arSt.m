%% plots the baseline results
dl_n = '../../Results2/CardiovascularSystem.mat'
dl_arst = '../../Results2/imp_arst.sdf'
% c:\Program Files\Dymola 2021x\bin\dsres2sdf.exe dl_n
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
dl_n = dymload(dl_n)
h5info(dl_arst)
%%
mmHg2SI = 133.322;
ml2SI = 1e-6;
bpm2SI = 1/60;
%%
time_n = dymget(dl_n, 'Time');
d = 1.6;
t_interval_n = [58.07, 58.2 + d]; % interval in seconds
td = t_interval_n(2) - t_interval_n(1)
i_int_n = (time_n >= t_interval_n(1) & time_n <= t_interval_n(2));
t_n = time_n - t_interval_n(1);

pb = dymget(dl_n, 'brachial_pressure')/mmHg2SI;
pbm = dymget(dl_n, 'brachial_pressure_mean')/mmHg2SI;
hr = dymget(dl_n, 'HR')/bpm2SI;
plv = dymget(dl_n, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
vlv = dymget(dl_n, 'V_LV')/ml2SI;
pow_lv_n = dymget(dl_n, 'heartComponent.ventricles.power_LV');


time_arst = h5read(dl_arst, '/Time');
t_interval_arst = [298.16, 298.16 + d]; % interval in seconds
i_int_arst = (time_arst >= t_interval_arst(1) & time_arst <= t_interval_arst(2));
t_arst = time_arst - t_interval_arst(1);

pbm_arst = h5read(dl_arst, '/brachial_pressure_mean')/mmHg2SI;
pb_arst = h5read(dl_arst, '/brachial_pressure')/mmHg2SI;
hr_arst = h5read(dl_arst, '/HR')/bpm2SI;
plv_arst = h5read(dl_arst, '/P_LV')/mmHg2SI;
vlv_arst = h5read(dl_arst, '/heartComponent/ventricles/V_LV')/ml2SI;
pow_lv_arst = h5read(dl_arst, '/heartComponent/ventricles/power_LV');

%
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 10);

% s_a1 = subplot(2, 2, 1);
hold on;
title('BP and HR of arterial stiffening');
plot(t_n(i_int_n), plv(i_int_n), 'g:', 'LineWidth', 1);
plot(t_arst(i_int_arst), plv_arst(i_int_arst), 'g', 'LineWidth', 0.5);


plot(t_n(i_int_n), pb(i_int_n), 'b:', 'LineWidth', 1);
plot(t_arst(i_int_arst), pb_arst(i_int_arst), 'b', 'LineWidth', 0.5);

plot(t_n(i_int_n), pbm(i_int_n), 'b:', 'LineWidth', 1.5);
plot(t_arst(i_int_arst), pbm_arst(i_int_arst), 'b--', 'LineWidth', 1.5);

plot(t_n(i_int_n), hr(i_int_n), 'r:', 'LineWidth', 1);
plot(t_arst(i_int_arst), hr_arst(i_int_arst), 'r', 'LineWidth', 0.5);


leg = legend('PLV N', 'PLV ArSt', 'BP N', 'BP ArSt', 'BPm N', 'BPm ArSt', 'HR N', 'HR ArSt', 'Location', 'SouthEast')
leg.ItemTokenSize = [15, 10];
xlabel('Time')
xlim([0, d]);

disp("HR change:" + num2str(mean(hr(i_int_n)) - mean(hr_arst(i_int_arst))) + " BPM");
% ylim([50, 120])

%% save figure
tw = 7;
th = 7;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th], 'PaperUnits', 'Centimeters', 'PaperSize', [tw, th])
%%
exportgraphics(gcf,'fig_R_ArSt.png','Resolution',300)
% print(gcf,'fig_R_base_300.png', '-dpng', '-r300')
% saveas(gcf,'fig1.svg')

%% PWV
sp = 5; % current speed in m/s

[a, vi_n] = max(plv(i_int_n));
[a, bi_n] = max(pb(i_int_n));

dtn = t_n(bi_n) - t_n(vi_n);
s = sp*dtn;


dist = 0.677*0.8;



% [a, vi_arst] = max(plv_arst(iarst_st:iarst_en));
% [a, bi_arst] = max(pb_arst(iarst_st:iarst_en))
% 
% dtarst = t_arst(bi_arst + iarst_st) - t_arst(vi_arst+iarst_st) ;
% sparst = s/dtarst
% 
pbc = h5read(dl_arst, '/carotid_pressure')/mmHg2SI;
pbf = h5read(dl_arst, '/femoral_pressure')/mmHg2SI;
time = h5read(dl_arst, '/Time');
% pbc = dymget(dl_n, 'carotid_pressure')/mmHg2SI;
% pbf = dymget(dl_n, 'femoral_pressure')/mmHg2SI;
% time = dymget(dl_n, 'Time');

iarst_en = size(time, 1);
iarst_st = find(time > time(iarst_en) - 1, 1);
inte = [iarst_st:iarst_en];

[a, vi_arst] = min(pbc(inte));
[a, bi_arst] = min(pbf(inte));
dtarst = time(bi_arst + iarst_st) - time(vi_arst+iarst_st) ;
sparst = dist/dtarst


clf;hold on;
plot(time(inte), pbc(inte), 'b');
plot(time(inte), pbf(inte), 'r')
plot(time(vi_arst+iarst_st), pbc(vi_arst+iarst_st), 'b*')
plot(time(bi_arst + iarst_st), pbf(bi_arst + iarst_st), 'r*')
title("PWV " + num2str(sparst) + "m/s")

% disp("Calculated PWV is " + num2str())
% plot(time_arst(inte), [diff(pbc(inte)); 0], '*');
% plot(time_arst(inte), [diff(pbf(inte)); 0], 'r*');
