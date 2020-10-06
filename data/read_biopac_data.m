
clear
close all;

%% BIG CITY lOOP
% files  = {'V_00_sup', 'V_01_sit', 'V_02_sit', 'V_03_sit', 'VEc_01_sup', 'VEc_02_sup', 'VEc_03_sup'};
files  = {'VEc_01_sup'};
folder = "Valsalva/";
% files = {'V_00_sit_01', 'V_00_sit_02', 'V_00_sit_03', 'V_00_sit_04'}
% files  = {'VEc_01_sup.mat'}
useSeqExt = true;
% folder = "";
% toSeconds = 1;
toSeconds = 60;
source_fs = 2000;
% source_fs = 200;
target_fs = 100;
figure(100);clf;hold on;
for i = 1:length(files)
    % filename = 'data/VEc_0' + string(i) + '_sup';
    % filename = 'Valsalva/V_0' + string(j) + '_sit';
    filename = char(files(i))
    % filename = "Valsalva/V_00_sup";
    filename_open = folder + filename + ".txt";
    % filename_open = folder + filename
    
    % dataset = import_txt_data('data/V_0' + string(i) + '_sit.txt');
    dataset = import_txt_data(filename_open);
    % open(filename_open)
    % dataset
    fs = 1/(dataset.Time(3)- dataset.Time(2))./toSeconds;
    
    
    %% decimate the dataset
    % decimation factor
%     d = 20;
    d = source_fs / target_fs;

    % smoothing factor
    s = 5*d;
    
    % fs = 2khz
    t = (dataset.Time)*toSeconds;
    bps = smooth( dataset.BP, s);
    tps = smooth(dataset.TP, s);
    
    % blood pressure smooth decimated
    bpsd = bps(1:d:end);
    % thoracic pressure smooth decimated
    tpsd = tps(1:d:end);
    
    % fs = 100hz
    td = t(1:d:end);
    fs = round(1/(td(2)-td(1)));
    
    %% find HR
    
    % peak cut-off
    pc = 5;
    %  get trend
    bpdst = smooth(bpsd, 100);
    % subtrack trend
    bpds = bpsd - bpdst-pc;
    bpds(bpds < 0 ) = 0;
    [~,maxs] = findpeaks(bpds,'MinPeakDistance',45);
    HR = [60;60./diff(td(maxs))];
    t_hr_beats = td(maxs);
    hr_all = interp1(td(maxs), HR, td);
    hr_all(isnan(hr_all)) = HR(1);
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
    plot(t_hr_beats, HR, '-*', 'linewidth', 2);
    plot(td, tpsd , 'linewidth', 2);
    
    title(filename, 'Interpreter', 'none')
    legend('Blood pressure', 'detected peaks', 'HR [bpm]', 'Thoracic pressure');
    
    %% find segments starts
    tp_treshhold = 20;
    % minimal distance of two consecutive peaks in seconds
    min_distance = 10;
    % minimal valsalva length
    min_length = 10;
    time_before = 20;
    time_after = 15 + 40;
    
    
    valsalva_on = tpsd > tp_treshhold;
    valsalva_starts = [0; valsalva_on(2:end) &  ~valsalva_on(1:end-1)];
    
    % rejecting starting points which do not continue after min_length
    valsalva_starts_positions = find(valsalva_starts == 1)
    vspl = length(valsalva_starts_positions)
    for j = 1:vspl
%         at time td(valsalva_starts_positions(j) + round(min_distance*fs))
        if tpsd(valsalva_starts_positions(j) + round(min_distance*fs)) < tp_treshhold
            valsalva_starts(valsalva_starts_positions(j)) = 0;
        end;
    end;
    
    % rejecting starting points too close to preceding one
    valsalva_starts_positions = find(valsalva_starts == 1)
    vspl = length(valsalva_starts_positions)
    for j = 2:vspl
        if td(valsalva_starts_positions(j-1)) + min_distance > td(valsalva_starts_positions(j))
            valsalva_starts(valsalva_starts_positions(j)) = 0;
        end;
    end;
    % now without rejected points
    valsalva_starts_positions = find(valsalva_starts == 1)
    
    plot(td(valsalva_starts_positions), tpsd(valsalva_starts_positions), 'b*', 'markersize', 10)
    plot(td(valsalva_starts_positions) - time_before, tpsd(valsalva_starts_positions), 'g*', 'markersize', 5)
    plot(td(valsalva_starts_positions) + time_after, tpsd(valsalva_starts_positions), 'r*', 'markersize', 5)
    %% cut sin segmentys and save
    for j = 1:length(valsalva_starts_positions)
        
        time_start = td(valsalva_starts_positions(j)) - time_before;
        time_end = time_start + time_after;
        incl = td > time_start & td < time_end;
        
        time = td(incl) - time_start;
        arterial_pressure = [time, bpsd(incl)];
        
        heart_rate = [time, hr_all(incl)];
        thoracic_pressure = [time, tpsd(incl)];
        if useSeqExt 
            segment = sprintf( '_%02d', j ) ;
        else
            segment = "";
        end
        save(folder + filename + segment + ".mat", 'arterial_pressure', 'heart_rate', 'thoracic_pressure')
        figure(100);hold on;
        plot(time, arterial_pressure(:, 2), time, heart_rate(:, 2), time, thoracic_pressure(:, 2))
        disp('Saved')
    end
end
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
