dataset = import_txt_data('Dan_lyingdown.txt');

% decimation factor
d = 20;

t = (dataset.Time - 7.09958)*60;
bp = dataset.BP*0.75;
bps = smooth(bp, 100);
% blood pressure smooth decimated
bpsd = bps(1:d:end);
td = t(1:d:end);

%%
%  get trend
    bpdst = smooth(bpsd, 100);
    % subtrack trend
    bpds = bpsd - bpdst-15;
    bpds(bpds < 0 ) = 0;
    
    
     [~,maxs] = findpeaks(bpds,'MinPeakDistance',55);
    
      HR = [60./diff(td(maxs));60];
     
     figure(2);
clf; hold on;
plot(t, bps, 'linewidth', 1);

% debug peak detection for HR
% plot(td, bpdst, 'linewidth', 2);
% plot(td, bpds, 'linewidth', 2);
plot(td(maxs), bpsd(maxs), 'rx');
plot(td(maxs), HR, 'linewidth', 2);

legend('Blood pressure, scaled by 0.75 [mmHg]', 'detected peaks', 'HR [bpm]');

%%     


%%
figure(1);
clf; hold on;
plot(td, bpsd, 'r-', 'linewidth', 1);
plot(t, bp, 'linewidth', 1);
tps = smooth(dataset.TP, 100);
plot(t, dataset.TP);
plot(t, tps, 'linewidth', 2);
%%
figure(2);
clf; hold on;
plot(t, [0; diff(dataset.BP)], 'linewidth', 2);




% plot(dataset.Time)
% filter out 20 Hz noise and smooth the TP data a bit

% staris with f of about 6 Hz

