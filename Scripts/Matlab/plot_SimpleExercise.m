clear;
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
run('../Figures/color_schema');

% Normal
[hrs_n,pas_n, pad_n, pcwp_n, co_n, pwrl_n, pwrr_n] = loadStepUpExercise('../../SEO_stepUp.mat');
% Identified to HFpEF
[hrs_hfpef,pas_hfpef, pad_hfpef, pcwp_hfpef, co_hfpef, pwrl_hfpef, pwrr_hfpef] = loadStepUpExercise('../../SEO_HFpEF_stepUp.mat');

%%

% data from our paper
% https://docs.google.com/spreadsheets/d/1mn0YtzfLdXnx4d4q5XqR5L8WkartYXXkaWIzZKXHbJo/edit#gid=926384656
hrs_d = [64; 153];
pas_nd = [120; 194];
pad_nd = [80; 81];
pcwp_nd = [8; 19.4];
% plvs_nd = [max(plv_n(int_1)); max(plv_n(int_2)); max(plv_n(int_3))]
co_nd = [5.87; 16.1];

% Data from Ot's patient nr7
hrs_hfd = [72; 96; 111];
pas_hfd = [153; 159; 163];
pad_hfd = [83; 82; 85];
pcwp_hfd = [13; 23; 25];
% plvs_nd = [max(plv_n(int_1)); max(plv_n(int_2)); max(plv_n(int_3))]
co_hfd = [5.76; 11.5; NaN];
%%
figure(1); clf;
subplot(221);hold on; cla;title('PA_{sys} and PA_{dia}');
plot(hrs_n, pas_n, '-*', 'Color', color_b);
plot(hrs_d, pas_nd, 'o', 'Color', color_b);
plot(hrs_n, pad_n, '--*', 'Color', color_b);
plot(hrs_d, pad_nd, 'o', 'Color', color_b);
plot(hrs_hfpef, pas_hfpef, '-*', 'Color', color_r);
plot(hrs_hfpef, pad_hfpef, '--*', 'Color', color_r);
plot(hrs_hfd, pas_hfd, 'o', 'Color', color_r);
plot(hrs_hfd, pad_hfd, 'o', 'Color', color_r);
xlabel('HR')

subplot(222);hold on; cla;title('PCWP mmHg');
plot(hrs_n, pcwp_n, '-*', 'Color', color_b);
plot(hrs_d, pcwp_nd, 'o', 'Color', color_b);
plot(hrs_hfpef, pcwp_hfpef, '-*', 'Color', color_r);
plot(hrs_hfd, pcwp_hfd, 'o', 'Color', color_r);
xlabel('HR')

subplot(223);cla;hold on;title('CO l/min');
plot(hrs_d, co_nd, 'o', 'Color', color_b);
plot(hrs_n, co_n, '-*', 'Color', color_b);
plot(hrs_hfd, co_hfd, 'o', 'Color', color_r);
plot(hrs_hfpef, co_hfpef, '-*', 'Color', color_r);
legend('Normal (data)', 'Normal (Model)', 'Patient nr7 (data)', 'Patient nr7 (model)','Location', 'best');
xlabel('HR')

subplot(224);hold on; cla;title('LV and RV power [W]');
plot(hrs_n, pwrl_n, '-*', 'Color', color_b);
plot(hrs_n, pwrr_n, '--*', 'Color', color_b);
plot(hrs_hfpef, pwrl_hfpef, '-*', 'Color', color_r);
plot(hrs_hfpef, pwrr_hfpef, '--*', 'Color', color_r);
xlabel('HR')

%%
T1 = table(hrs_n,pas_n, pad_n, pcwp_n, co_n, pwrl_n, pwrr_n)
T2 = table(hrs_hfpef,pas_hfpef, pad_hfpef, pcwp_hfpef, co_hfpef, pwrl_hfpef, pwrr_hfpef)
T3 = table(hrs_d,pas_nd, pad_nd, pcwp_nd,co_nd)
T4 = table(hrs_hfd,pas_hfd, pad_hfd, pcwp_hfd,co_hfd)
writetable(T3,filename,'Sheet','Patient (data)')

filename = 'patientdata.xlsx';
writetable(T1,filename,'Sheet','Normal {model)')
writetable(T3,filename,'Sheet','Normal (data)')
writetable(T2,filename,'Sheet','Patient (model)')
writetable(T4,filename,'Sheet','Patient (data)')
%%
function [hrs,pas, pad, pcwp, co, pwrl, pwrr] = loadStepUpExercise(datafile)
    run('../Figures/color_schema');
    dl_n = dymload(datafile);

    steps_up = dymget(dl_n, 'exercise.y_');
    t_n = dymget(dl_n, 'Time');
    hr_n = dymget(dl_n, 'HR');
    pa = dymget(dl_n, 'brachial_pressure')/mmHg2SI;
    pwr_lv = dymget(dl_n, 'heartComponent.ventricles.power_LV');
    pwr_rv = dymget(dl_n, 'heartComponent.ventricles.power_RV');
    pcwp = dymget(dl_n, 'P_pv')/mmHg2SI;
%     plv = dymget(dl_n, 'heartComponent.ventricles.P_LV')/mmHg2SI;
    co = dymget(dl_n, 'CO')*1000*60;

    % get times
    times_i = find(diff(steps_up) > 0);
    total = size(steps_up, 1) -1;
    times_i = [times_i; total];
    hrs = hr_n(times_i)*60;

    % 3 secs of moving averaging window
    win_len = find(t_n > 3, 1);
    pas = movmax(pa, [win_len, 0]);
    pad = movmin(pa, [win_len, 0]);
    pcwp = movmean(pcwp, [win_len, 0]);
    co = movmean(co, [win_len, 0]);
    pwrl = movmax(pwr_lv, [win_len, 0]);
    pwrr = movmax(pwr_rv, [win_len, 0]);
    
    pas = pas(times_i);
    pad = pad(times_i);
    pcwp = pcwp(times_i);
    co = co(times_i);
    pwrl = pwrl(times_i);
    pwrr = pwrr(times_i);
end

% function avg = movavg(data, window)
%     win_len = window(2) - window(1);
%     for i = window(1):size(data, 1),
%         data(i-window(1) + 1:window)
        