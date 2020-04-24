% plots all the values at once to have avg reaction
clear
close all
% sitting
% files = {"V_00_sit_01", "V_00_sit_02", "V_00_sit_03", "V_00_sit_04", "V_01_sit_01", "V_01_sit_02", "V_01_sit_03", "V_02_sit_01", "V_02_sit_02", "V_03_sit_01", "V_03_sit_02"}
% supine
% all
% files = {"V_00_sup", "VEc_01_sup_01", "VEc_01_sup_02", "VEc_01_sup_03", "VEc_02_sup_01", "VEc_02_sup_02", "VEc_03_sup_01"}
% subj
files = {"VEc_01_sup_01", "VEc_01_sup_02", "VEc_01_sup_03"}

folder = "Valsalva/";
%% plot it
figure(100);
clf; hold on;
tp_summa = []
bp_summa = []
hr_summa = []
fs = 100
for file = files
    % vars: 'time', 'arterial_pressure', 'heart_rate', 'thoracic_pressure'
    load(folder + string(file) + ".mat")
    time = arterial_pressure(:, 1);
    bp_mean = smooth(arterial_pressure(:, 2), 1*fs);
%     clf;plot(time, arterial_pressure, time, bp_mean, 'linewidth', 2)
    
    tp_summa = AddTOSum(tp_summa, thoracic_pressure(:, 2) );
    bp_summa = AddTOSum(bp_summa, bp_mean);
    hr_summa = AddTOSum(hr_summa, heart_rate(:, 2) );
%     add = arterial_pressure
%     
%     if length(base) > length(base)
%         base2 = base
%         add2 = add
%     else
%         base2 = add
%         add2 = base
figure();hold on;
    plot(time, arterial_pressure(:, 2), 'linewidth', 0.5);
    plot(time, bp_mean, 'linewidth', 2);
    plot(time, heart_rate(:, 2), 'linewidth', 2);
    plot(time, thoracic_pressure(:, 2) , 'linewidth', 3);
    title(file, 'Interpreter', 'none')
    legend('Blood pressure', 'HR', 'Thoracic pressure');

    figure(100);
    plot(time, bp_mean, 'linewidth', 1);
     plot(time, arterial_pressure(:, 2), 'linewidth', 1);
    plot(time, heart_rate(:, 2), 'linewidth', 1);
    plot(time, thoracic_pressure(:, 2) , 'linewidth', 1);
end

    tp_mean = tp_summa(1:length(time)) / length(files);
    bp_mean = bp_summa(1:length(time)) / length(files);
    hr_mean = hr_summa(1:length(time)) / length(files);

    plot(time, tp_mean, 'linewidth', 4);
    plot(time, bp_mean , 'linewidth', 4);
    plot(time, hr_mean, 'linewidth', 4);
    