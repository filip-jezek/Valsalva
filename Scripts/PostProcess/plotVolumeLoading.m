
close all;
clear;

% (D)ata (P)reserved ejection fraction by adjusting the (k) stiffness slope
Dpk = readtable('VolumeLoading_HFpEF_baro.csv',  'PreserveVariableNames',false);
% (D)ata (P)reserved ejection fraction by (d)dilating the heart's AmRef
Drd = readtable('VolumeLoading_HFpEFDilated_baro.csv',  'PreserveVariableNames',false);
% (D)ata (r)reduced ejection fraction by adjusting the LV (c)ontractility
Drc = readtable('VolumeLoading_CHF_baro.csv',  'PreserveVariableNames',false);


% Dpk.Properties.VariableNames
%%
figure(1);
tiledlayout('flow');clf;

% Volume
% figure();
% nexttile;
% hold on;xlabel('PCWP mmHg');
% title('volume to compensation (Liters)')
% plot(Dpk.PCWP_mmHg_, Dpk.vol_L_, 'b*')
% plot(Drd.PCWP_mmHg_, Drd.vol_L_, 'k*')
% plot(Drc.PCWP_mmHg_, Drc.vol_L_, 'r*')
% legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')
nexttile;
hold on;xlabel('PCWP mmHg');
title('Excess volume in IS (Liters)')
plot(Dpk.PCWP_mmHg_, Dpk.ExcessInterstitialVolume_L_, 'b*')
plot(Drd.PCWP_mmHg_, Drd.ExcessInterstitialVolume_L_, 'k*')
plot(Drc.PCWP_mmHg_, Drc.ExcessInterstitialVolume_L_, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% CO
nexttile;hold on;xlabel('PCWP mmHg');
title('CO lpm')
plot(Dpk.PCWP_mmHg_, Dpk.CO_L_min_, 'b*')
plot(Drd.PCWP_mmHg_, Drd.CO_L_min_, 'k*')
plot(Drc.PCWP_mmHg_, Drc.CO_L_min_, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% EF
nexttile;hold on;xlabel('PCWP mmHg');
title('EF')
plot(Dpk.PCWP_mmHg_, Dpk.EF, 'b*')
plot(Drd.PCWP_mmHg_, Drd.EF, 'k*')
plot(Drc.PCWP_mmHg_, Drc.EF, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')


% eGFR x PCWP
nexttile;hold on;xlabel('PCWP mmHg');
title('eGFR')
plot(Dpk.PCWP_mmHg_, Dpk.eGFR, 'b*')
plot(Drd.PCWP_mmHg_, Drd.eGFR, 'k*')
plot(Drc.PCWP_mmHg_, Drc.eGFR, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% HR
nexttile;hold on;xlabel('PCWP mmHg');
title('HR')
plot(Dpk.PCWP_mmHg_, Dpk.HR, 'b*')
plot(Drd.PCWP_mmHg_, Drd.HR, 'k*')
plot(Drc.PCWP_mmHg_, Drc.HR, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% SV
nexttile;hold on;xlabel('PCWP mmHg');
title('SV')
plot(Dpk.PCWP_mmHg_, Dpk.EF.*Dpk.EDVMl./100 , 'b*')
plot(Drd.PCWP_mmHg_, Drd.EF.*Drd.EDVMl./100 , 'k*')
plot(Drc.PCWP_mmHg_, Drc.EF.*Drc.EDVMl./100 , 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

%%
figure(2);
tiledlayout('flow');clf;

% Volume
% figure();
% nexttile;
% hold on;xlabel('Central venous pressure mmHg');
% title('volume to compensation (Liters)')
% plot(Dpk.PCWP_mmHg_, Dpk.vol_L_, 'b*')
% plot(Drd.PCWP_mmHg_, Drd.vol_L_, 'k*')
% plot(Drc.PCWP_mmHg_, Drc.vol_L_, 'r*')
% legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')
nexttile;
hold on;xlabel('Central venous pressure mmHg');
title('Excess volume in IS (Liters)')
plot(Dpk.CVP, Dpk.ExcessInterstitialVolume_L_, 'b*')
plot(Drd.CVP, Drd.ExcessInterstitialVolume_L_, 'k*')
plot(Drc.CVP, Drc.ExcessInterstitialVolume_L_, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% CO
nexttile;hold on;xlabel('Central venous pressure mmHg');
title('CO lpm')
plot(Dpk.CVP, Dpk.CO_L_min_, 'b*')
plot(Drd.CVP, Drd.CO_L_min_, 'k*')
plot(Drc.CVP, Drc.CO_L_min_, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% EF
nexttile;hold on;xlabel('Central venous pressure mmHg');
title('EF')
plot(Dpk.CVP, Dpk.EF, 'b*')
plot(Drd.CVP, Drd.EF, 'k*')
plot(Drc.CVP, Drc.EF, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')


% eGFR x PCWP
nexttile;hold on;xlabel('Central venous pressure mmHg');
title('eGFR')
plot(Dpk.CVP, Dpk.eGFR, 'b*')
plot(Drd.CVP, Drd.eGFR, 'k*')
plot(Drc.CVP, Drc.eGFR, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% HR
nexttile;hold on;xlabel('Central venous pressure mmHg');
title('HR')
plot(Dpk.CVP, Dpk.HR, 'b*')
plot(Drd.CVP, Drd.HR, 'k*')
plot(Drc.CVP, Drc.HR, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

% SV
nexttile;hold on;xlabel('Central venous pressure mmHg');
title('SV')
plot(Dpk.CVP, Dpk.EF.*Dpk.EDVMl./100 , 'b*')
plot(Drd.CVP, Drd.EF.*Drd.EDVMl./100 , 'k*')
plot(Drc.CVP, Drc.EF.*Drc.EDVMl./100 , 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

%%

figure(3);
tiledlayout('flow');clf;

% Volume
% figure();
nexttile;
hold on;xlabel('Central venous pressure mmHg');
title('volume to compensation (Liters)')
plot(Dpk.f_LV, Dpk.vol_L_, 'b*')
plot(Drd.f_LV, Drd.vol_L_, 'k*')
plot(Drc.f_LV, Drc.vol_L_, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

nexttile; hold on;title('Excess volume in IS (Liters)');xlabel('Volume loaded (L)');
plot(Dpk.vol_L_, Dpk.ExcessInterstitialVolume_L_, 'b*')
plot(Drd.vol_L_, Drd.ExcessInterstitialVolume_L_, 'k*')
plot(Drc.vol_L_, Drc.ExcessInterstitialVolume_L_, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

nexttile; hold on;title('CVP on PCWP (L)');xlabel('PCWP (Liters)'); 
plot(Dpk.CVP, Dpk.PCWP_mmHg_, 'b*')
plot(Drd.CVP, Drd.PCWP_mmHg_, 'k*')
plot(Drc.CVP, Drc.PCWP_mmHg_, 'r*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

nexttile; hold on;title('ISF on f\_LV');xlabel('f_lv -'); 
plot(Dpk.f_LV, Dpk.ExcessInterstitialVolume_L_, 'b-*')
plot(Drd.f_LV, Drd.ExcessInterstitialVolume_L_, 'k-*')
plot(Drc.f_LV, Drc.ExcessInterstitialVolume_L_, 'r-*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

nexttile; hold on;title('Time on f\_LV');xlabel('f_lv -'); 
plot(Dpk.f_LV, Dpk.time_s_, 'b-*')
plot(Drd.f_LV, Drd.time_s_, 'k-*')
plot(Drc.f_LV, Drc.time_s_, 'r-*')
legend('HFpEF, stiffness', 'HFrEF, dilated', 'HFrEF, contractility')

