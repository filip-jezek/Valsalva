% work in progrezs
function s = plotExSubplot(filename, tit, fakeIt, drawLegend)
if fakeIt
    title(tit)
    ylabel('zadky')
    xlabel('P5edky')
    
    yyaxis right;
    ylabel('spodky')
    s.pow_lv = rand(11,1);
    return
end;
    
%%
% filename = '../../Results2/ExStep_base.mat'
% figure(3); clf;hold on;
%%
dl = dymload(char(filename));
color_schema;

time = dymget(dl, 'Time');
% t_interval_n = [58.13, 58.13 + 1.6]; % interval in seconds
% td_n = t_interval_n(2) - t_interval_n(1)
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
% i_int_n = (time >= t_interval_n(1) & time <= t_interval_n(2));

% t_n = time - t_interval_n(1);
pa_m = dymget(dl, 'brachial_pressure_mean')/mmHg2SI;
pa_s = dymget(dl, 'brachial_pressure_systolic')/mmHg2SI;
pa_d = dymget(dl, 'brachial_pressure_diastolic')/mmHg2SI;
pw = dymget(dl, 'P_pv')/mmHg2SI;
pvc = dymget(dl, 'P_sv')/mmHg2SI;

el = dymget(dl, 'Exercise.y_')*100;

pow_lv = dymget(dl, 'heartComponent.ventricles.power_LV');
% pow_rv = h5read(dl, '/heartComponent/ventricles/power_RV');
vlv = dymget(dl, 'V_LV')/ml2SI;

co = dymget(dl, 'CO')/mlPmin2SI;
q_ex = dymget(dl, 'q_exercised_avg')/mlPmin2SI;

try
    % try corvia data
    q_corv = - dymget(dl, 'r_SystemicVenousInflow.volumeFlowRate')/mlPmin2SI;
catch KE
    q_corv = zeros(size(time));
end
%%
ts = 35;
te = 36.6;
td = te - ts;
diffs = [1; find( diff(el) > 0)];

set_ef = [];
set_edv = [];
set_pow_lvs = [];
set_pa_m = [];
set_pa_s = [];
set_pa_d = [];
set_vlv = [];
set_pw = [];
set_pvc = [];
set_q_corv = [];
% set_pow_rv = [];
set_co = [];
set_q_ex = [];


for i = 1:size(diffs)
%      i = 9

    i_s = find(time > time(diffs(i)) + ts, 1);
    i_e = find(time > time(diffs(i)) + te, 1);
%     t = time - time(i_s);
    i_int = [i_s:i_e];
    if i == 9
        t = time - time(i_s);
        i_int90 = i_int;
    end
%     plot(vlv(i_s:i_e), plv(i_s:i_e))
%     plot(vrv(i_s:i_e), prv(i_s:i_e))
%     plot(vla(i_s:i_e), pla(i_s:i_e))
    
    sv =  max(vlv(i_s:i_e)) -  min(vlv(i_s:i_e));
    set_ef = [set_ef, sv / max(vlv(i_s:i_e))];
    set_edv = [set_edv, min(vlv(i_s:i_e))];
    set_pa_m = [set_pa_m, mean(pa_m(i_int))];
    set_pa_s = [set_pa_s, mean(pa_s(i_int))];
    set_pa_d = [set_pa_d, mean(pa_d(i_int))];
    set_pw = [set_pw, mean(pw(i_int))];
    set_pow_lvs = [set_pow_lvs , pow_lv(i_s)];
%     pow_rvs = [pow_rvs , pow_rv(i_s)];
    set_co = [set_co, mean(co(i_int))];
    set_q_ex = [set_q_ex, mean(q_ex(i_int))];
    set_q_corv = [set_q_corv, mean(q_corv(i_int))];
    set_pvc = [set_pvc, mean(pvc(i_int))];
end

s.ef = set_ef;
s.pow_lv = set_pow_lvs;
s.q_c = set_q_corv;
s.pvc = set_pvc;
s.pcwp = set_pw;
s.co = set_co;
s.pa = set_pa_m
% plot([0:10:100], edv, 'b*-')
% plot([0:10:100], pow_rvs, 'r*-')
% plot([0:10:100], pow_lvs, 'b*-')
% plot([0:10:100], ef*10, 'b*-')
% plot([0:10:100], set_co, 'b*-')
% plot([0:10:100], set_q_ex, 'r*-')
%%
hold on;
title(tit, 'interpreter', 'None')
t_ax = [0:10:100];
gf = fill([t_ax, fliplr(t_ax)], [set_pa_s, fliplr(set_pa_d)], color_lb,'EdgeColor',color_lb);

gpm = plot(t_ax, set_pa_m, 'x-', 'LineWidth', 1, 'Color', color_b);
% gpw = plot(t_ax, set_pw*10, '--', 'Color', color_g, 'LineWidth', 1);
% plot(t_ax, set_pa_s, '.', 'LineStyle', 'None', 'Color', color_b);
% plot(t_ax, set_pa_d, '.', 'LineStyle', 'None', 'Color', color_b);
ylim([0, 300])
ylabel('Pressure (mmHg)')

yyaxis right;
gco = plot(t_ax, set_q_ex, 'x:', 'LineWidth', 1, 'Color', color_m);
gqe = plot(t_ax, set_co, 'x-', 'LineWidth', 1, 'Color', color_m);
ylim([0, 30]);
g = gca; g.YAxis(2).Color = color_m;
g.TickLength = [0.05, 0.05];

xlabel('Exercise (% of max)');
ylabel('Blood flow (L/min)');
if drawLegend
%     leg = legend([gf, gpm, gco, gqe, gpw], 'PA', 'PA mean', 'CO', 'Q Ex','PPW (x10)', 'Location', 'NorthEast', 'Orientation', 'Horizontal');
    leg = legend([gf, gpm, gco, gqe], 'PA', 'PA mean', 'CO', 'Q Ex','Location', 'NorthEast', 'Orientation', 'Horizontal');
    leg.ItemTokenSize = [20, 150];
end
s.t_ax = t_ax;