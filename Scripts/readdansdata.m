% dataset = import_txt_data('data/Dan_lyingdown.txt');
% dataset = import_txt_data('data/Filip sitting.txt');
% dataset = import_txt_data('data/Noah sitting.txt');
% dataset = import_txt_data('data/Filip with ECHO.txt');
clear
i = 3
% filename = 'data/VEc_0' + string(i) + '_sup';
filename = 'data/V_0' + string(i) + '_sit';
filename_open = filename + '.txt'
% dataset = import_txt_data('data/V_0' + string(i) + '_sit.txt');
dataset = import_txt_data(filename_open);


%% decimate the dataset
% decimation factor
d = 20;
% smoothing factor
s = 100;

% fs = 2khz
t = (dataset.Time)*60;
bps = smooth( dataset.BP, 100);
tps = smooth(dataset.TP, 100);

% blood pressure smooth decimated
bpsd = bps(1:d:end);
% thoracic pressure smooth decimated
tpsd = tps(1:d:end);

% fs = 100hz
td = t(1:d:end);

%% find HR

% peak cut-off
pc = 10;
%  get trend
bpdst = smooth(bpsd, s);
% subtrack trend
bpds = bpsd - bpdst-pc;
bpds(bpds < 0 ) = 0;
[~,maxs] = findpeaks(bpds,'MinPeakDistance',45);
HR = [60./diff(td(maxs));60];
%%     
figure(10 + i);
clf; hold on;
% plot(t, bp, 'linewidth', 1);
% plot(t, bp, 'linewidth', 1);
plot(td, bpsd, 'linewidth', 1);

% % debug peak detection for HR
plot(td, bpdst, 'linewidth', 2);
plot(td, bpds, 'linewidth', 2);

plot(td(maxs), bpsd(maxs), 'rx');
plot(td(maxs), HR, 'linewidth', 2);
plot(td, tpsd , 'linewidth', 2);

legend('Blood pressure', 'detected peaks', 'HR [bpm]', 'Thoracic pressure');

%% save
time_start = 0;
time_end = td(end);
incl = td > time_start & td < time_end;
time = td(incl);
bp = bpsd(incl);
hr_all = interp1(td(maxs), HR, td);
hr = hr_all(incl);
tp = tpsd(incl);
save(filename + '.mat', 'time', 'bp', 'hr', 'tp')
disp('Saved')

% %%
% figure(2);
% clf; hold on;
% plot(td, bpsd, 'r-', 'linewidth', 1);
% plot(t, bp, 'linewidth', 1);
% tps = smooth(dataset.TP, 100);
% plot(t, dataset.TP);
% plot(t, tps, 'linewidth', 2);
% %%
% figure(2);
% clf; hold on;
% plot(t, [0; diff(dataset.BP)], 'linewidth', 2);
% 
% 
% 
% 
% % plot(dataset.Time)
% % filter out 20 Hz noise and smooth the TP data a bit
% 
% % staris with f of about 6 Hz
% 
