%% read it from excel
bpd = readallsheets('Valsalva/BP_data.xlsx');
tpd = readallsheets('Valsalva/Valsalva_data');
hrd = readallsheets('Valsalva/HR_data');

%%
% k = 2;

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
% figure(1);clf;hold on;
% plot(tpd{k}(:, 1), tpd{k}(:, 2), 'r');
% plot(tpd{k}(:, 1), tps, 'b', 'LineWidth', 2);

%% produce the output
    tps = smooth(tpd{k}(:, 2), 20);
    thoracic_pressure = [tpd{k}(:, 1), tps];
    arterial_pressure = [bpd{k}(:, 1), bpd{k}(:, 2)];
    heart_rate = [hrd{k}(:, 1), hrd{k}(:, 5)];
    save(['valsalva_experiment', num2str(k),  '.mat'], 'thoracic_pressure', 'arterial_pressure', 'heart_rate');

%% save as txt for further processing
% saving in minutes to be consistent
time = bpd{k}(:, 1);
tp = interp1(tpd{k}(:, 1), tps, time, 'nearest');
tp(isnan(tp))=tps(end);
bp = bpd{k}(:, 2);

dataset = [time';bp';tp'];
% example line> 0.00645833	14.6472	1.97803	
% sprintf("%8.8f\t%4.4f\t%2.5f", a)
file = fopen("Valsalva/V_00_sit_0" + string(k) + ".txt",'w');
fprintf(file, "%8.8f	%4.4f	%2.5f	\r\n", dataset);
fclose(file)

end;    