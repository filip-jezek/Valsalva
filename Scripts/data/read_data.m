%% read it from excel
bpd = readallsheets('data/BP_data.xlsx');
tpd = readallsheets('data/Valsalva_data');
hrd = readallsheets('data/HR_data');

%%
k = 2;

%% the loop
for k = 1:min([numel(tpd), numel(bpd), numel(hrd)])
%% find BP peaks
    

    bpdsm = smooth(bpd{k}(:, 1), bpd{k}(:, 2), 0.01);
    bpds = bpd{k}(:, 2) - bpdsm;
    bpds(bpds < 0 ) = 0;
    [~,maxs] = findpeaks(bpds,'MinPeakDistance',55);
    
     HR = [60./diff(bpd{k}(maxs, 1));0];
%  test plot the HR

    

    
    figure(k);clf;hold on;    
%     set(gca,'ColorOrderIndex',k);
    plot(bpd{k}(:, 1), bpd{k}(:, 2), '-');
    plot(bpd{k}(:, 1), bpds, 'k-');
    plot(bpd{k}(:, 1), bpdsm, 'k-');
    plot(bpd{k}(maxs, 1), bpd{k}(maxs, 2), 'ro', 'MarkerSize', 12);
%%    
    set(gca,'ColorOrderIndex',k);
    plot(hrd{k}(:, 1), hrd{k}(:, 2), '*', 'MarkerSize', 12)
    set(gca,'ColorOrderIndex',k);
    plot(hrd{k}(:, 1), hrd{k}(:, 3), '+', 'MarkerSize', 12)
    plot(hrd{k}(:, 1), hrd{k}(:, 5), 'r-', 'LineWidth', 1);
    plot(bpd{k}(maxs, 1)+0.3, HR, 'g-', 'LineWidth', 1);
%     tps = smooth(tpd{k}(:, 1), tpd{k}(:, 2), 0.02, 'loess');
    tps = smooth(tpd{k}(:, 2), 20);
    plot(tpd{k}(:, 1), tps, 'k-', 'LineWidth', 2);


%% smooth the thoracic pressure
% tps = smooth(tpd{k}(:, 1), tpd{k}(:, 2), 0.04, 'loess');
tps = smooth(tpd{k}(:, 2), 20);
figure(1);clf;hold on;
plot(tpd{k}(:, 1), tpd{k}(:, 2), 'r');
plot(tpd{k}(:, 1), tps, 'b', 'LineWidth', 2);

%% produce the output
    tps = smooth(tpd{k}(:, 2), 20);
    thoracic_pressure = [tpd{k}(:, 1), tps];
    arterial_pressure = [bpd{k}(:, 1), bpd{k}(:, 2)];
    heart_rate = [hrd{k}(:, 1), hrd{k}(:, 5)];
    save(['valsalva_experiment', num2str(k),  '.mat'], 'thoracic_pressure', 'arterial_pressure', 'heart_rate');
end;    