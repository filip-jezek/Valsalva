
% (D)ata (P)reserved ejection fraction by adjusting the (k) stiffness slope
Dpk = readtable('VolumeLoading_HFpEF_baro.csv',  'PreserveVariableNames',false);
% (D)ata (P)reserved ejection fraction by (d)dilating the heart's AmRef
Drd = readtable('VolumeLoading_HFpEFDilated_baro.csv',  'PreserveVariableNames',false);
% (D)ata (r)reduced ejection fraction by adjusting the LV (c)ontractility
Drc = readtable('VolumeLoading_CHF_baro.csv',  'PreserveVariableNames',false);


% Dpk.Properties.VariableNames

close all;

tiledlayout('flow');clf;

% Volume
% figure();
nexttile;
hold on;xlabel('PCWP mmHg');
title('volume to compensation (Liters)')
plot(Dpk.PCWP_mmHg_, Dpk.vol_L_, 'b*')
plot(Drd.PCWP_mmHg_, Drd.vol_L_, 'k*')
plot(Drc.PCWP_mmHg_, Drc.vol_L_, 'r*')
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