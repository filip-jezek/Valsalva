% loads and Plots live starling curve volume loading result
close all;
clear;
% set_norm = importVolLoadResultFile("..\..\Results\baro_VolLoad_result_normal_.csv");
% set_norm.glyph = '*';
% set_hfr = importVolLoadResultFile("..\..\Results\baro_VolLoad_result_hfref_.csv");
% set_hfp = importVolLoadResultFile("..\..\Results\baro_VolLoad_result_hfpef_.csv");
% set_norm_nb = importVolLoadResultFile("..\..\Results\nobaro_VolLoad_result_normal_.csv");
% set_hfr_nb = importVolLoadResultFile("..\..\Results\nobaro_VolLoad_result_hfref_.csv");
% set_hfp_nb = importVolLoadResultFile("..\..\Results\nobaro_VolLoad_result_hfpef_.csv");

sethfr25.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_nobaro_LV25_.csv");
sethfr25.glyph = '.:';
sethfr25.title = 'HFr25 NB';
sethfr50.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_nobaro_LV50_.csv");
sethfr50.glyph = '.:';
sethfr50.title = 'HFr50 NB';
setnorm.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_nobaro_LV100_.csv");
setnorm.glyph = '.:';
setnorm.title = 'Norm NB';

sethfr25b.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_baro_LV25_.csv");
sethfr25b.glyph = 'd-';
sethfr25b.title = 'HFr25 B';
sethfr50b.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_baro_LV50_.csv");
sethfr50b.glyph = 's-';
sethfr50b.title = 'HFr50 B';
setnormb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_baro_LV100_.csv");
setnormb.glyph = 'o-';
setnormb.title = 'Norm B';


setnorml.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_nobaro_Lumens_LV100_.csv");
setnorml.glyph = '*:';
setnorml.title = 'Lumens NB';
sethfr50l.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_nobaro_Lumens_LV50_.csv");
sethfr50l.glyph = '+:';
sethfr50l.title = 'Lumens HFr50 NB';
sethfr25l.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_nobaro_Lumens_LV25_.csv");
sethfr25l.glyph = '+:';
sethfr25l.title = 'Lumens HFr25 NB';

setnormadj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSadj_nobaro_LV100_.csv");
setnormadj.glyph = '*:';
setnormadj.title = 'Adjusted NB';
sethfr50adj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSadj_nobaro_LV50_.csv");
sethfr50adj.glyph = '+:';
sethfr50adj.title = 'Adjusted HFr50 NB';
sethfr25adj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSadj_nobaro_LV25_.csv");
sethfr25adj.glyph = '+:';
sethfr25adj.title = 'Adjusted HFr25 NB';

setnormadjb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSadj_baro_LV100_.csv");
setnormadjb.glyph = 'o-';
setnormadjb.title = 'Adjusted B';
sethfr50adjb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSadj_baro_LV50_.csv");
sethfr50adjb.glyph = 'o-';
sethfr50adjb.title = 'Adjusted HFr50 B';
sethfr25adjb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSadj_baro_LV25_.csv");
sethfr25adjb.glyph = 'o-';
sethfr25adjb.title = 'Adjusted HFr25 B';

setnormSLradj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_nobaro100_.csv");
setnormSLradj.glyph = '*:';
setnormSLradj.title = 'SLr Adjusted NB';
sethfr50SLradj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_nobaro50_.csv");
sethfr50SLradj.glyph = '+:';
sethfr50SLradj.title = 'SLr Adjusted HFr50 NB';
sethfr30SLradj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_nobaro30_.csv");
sethfr30SLradj.glyph = '+:';
sethfr30SLradj.title = 'SLr Adjusted HFr30 NB';
sethfrR50SLradj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_nobaro_RV50_.csv");
sethfrR50SLradj.glyph = '+--';
sethfrR50SLradj.title = 'SLr Adjusted HFr RV50 NB';
sethfrR30SLradj.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_nobaro_RV30_.csv");
sethfrR30SLradj.glyph = '+--';
sethfrR30SLradj.title = 'SLr Adjusted HFr RV30 NB';

setnormSLradjb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_baro_LV100_.csv");
setnormSLradjb.glyph = 'o-';
setnormSLradjb.title = 'SLr Adjusted B';
sethfr50SLradjb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_baro_LV50_.csv");
sethfr50SLradjb.glyph = 'o-';
sethfr50SLradjb.title = 'SLr Adjusted HFr50 B';
sethfr30SLradjb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVSslr_baro_LV30_.csv");
sethfr30SLradjb.glyph = 'o-';
sethfr30SLradjb.title = 'SLr Adjusted HFr30 B';


setnormSLrOpt.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_nobaroSlrOpt_LV100_.csv");
setnormSLrOpt.glyph = 'x:';
setnormSLrOpt.title = 'SLr optimized NB';


setnormSLrOptb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_baroSlrOpt_LV100_.csv");
setnormSLrOptb.glyph = 'o-';
setnormSLrOptb.title = 'SLr optimized B';

setnormSLrOptnb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_baroSlrOptNatBaro_LV100_.csv");
setnormSLrOptnb.glyph = 'v-';
setnormSLrOptnb.title = 'SLr optimized naturalB';

setOlsonsnb.tab = importVolLoadResultFile("..\..\Results\VolLoad_result_CVS_Ols_nobaro_LV100_.csv");
setOlsonsnb.glyph = '*:';
setOlsonsnb.title = 'Olsons NB';

%% Compensated LV failure
setCompensatingLVc.tab = importVolLoadResultFile("..\..\Results\VolumeLoading_LV.csv");
setCompensatingLVc.glyph = '*-';
setCompensatingLVc.title = 'HC';
var_set = {setCompensatingLVc};
prenom = '';

%% Exercise at compensated HF

setLV10CompensatedExe.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_LV10.csv");
setLV10CompensatedExe.glyph = '*-';
setLV10CompensatedExe.title = '0.1 C_{LV}+EXE';

setLV20CompensatedExe.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_LV20.csv");
setLV20CompensatedExe.glyph = '*-';
setLV20CompensatedExe.title = '0.2 C_{LV}+EXE';

setLV40CompensatedExe.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_LV40.csv");
setLV40CompensatedExe.glyph = '*-';
setLV40CompensatedExe.title = '0.4 C_{LV}+EXE';

setLV60CompensatedExe.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_LV60.csv");
setLV60CompensatedExe.glyph = '*-';
setLV60CompensatedExe.title = '0.6 C_{LV}+EXE';

setLV80CompensatedExe.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_LV80.csv");
setLV80CompensatedExe.glyph = '*-';
setLV80CompensatedExe.title = '0.8 C_{LV}+EXE';

setLV100CompensatedExe.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_LV100.csv");
setLV100CompensatedExe.glyph = '*-';
setLV100CompensatedExe.title = '1.0 C_{LV}+EXE';

var_set = {setLV10CompensatedExe, setLV20CompensatedExe, setLV40CompensatedExe, setLV60CompensatedExe, setLV80CompensatedExe, setLV100CompensatedExe};
var_set = {setLV10CompensatedExe, setLV20CompensatedExe, setLV60CompensatedExe, setLV100CompensatedExe};
checkAllVarSets(var_set)
%% Compansated Exercise by failure

setCompensatedExe0.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_Exe0.mat.csv");
setCompensatedExe0.glyph = '*-';
setCompensatedExe0.title = 'LV_{sys, fail} Exe 0';

setCompensatedExe40.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_Exe40.mat.csv");
setCompensatedExe40.glyph = '*-';
setCompensatedExe40.title = 'LV_{sys, fail} Exe 40';

setCompensatedExe80.tab = importVolLoadResultFile("..\..\Results\VolLoad_CVS_baroExe_Exe80.mat.csv");
setCompensatedExe80.glyph = '*-';
setCompensatedExe80.title = 'LV_{sys, fail} Exe 80';
var_set = {setCompensatedExe0, setCompensatedExe40, setCompensatedExe80};
%% Setting var_set variants
var_set = {sethfr25, sethfr25b,sethfr50, sethfr50b,setnorm, setnormb};

var_set = {setnormb, sethfr50b,sethfr25b};

var_set = {setnorm,setnorml, setnormSLradj};

var_set = {setnorm,setnormb,setnormadj, setnormadjb};
var_set = {setnorm,setnormb};
var_set = {setnormadjb, sethfr50adjb, sethfr25adjb};
var_set = {setnormadj, sethfr50adj, sethfr25adj};
var_set = {setnormSLradj, sethfr50SLradj, sethfr30SLradj};
var_set = {setnormSLradj, setnormSLradjb,sethfr30SLradj, sethfr30SLradjb};
var_set = {setnormSLradj, setnormSLradjb,sethfr50SLradj, sethfr50SLradjb, sethfr30SLradj, sethfr30SLradjb};
var_set = {setnormSLradj, sethfr50SLradj, sethfrR50SLradj, sethfr30SLradj, sethfrR30SLradj};
% var_set = {setnormadjb,sethfr25adjb, sethfr25adj};
var_set = {sethfr25b, sethfr50b, setnormb, setnormadjb, sethfr50adjb, sethfr25adjb, setnormSLradjb, sethfr50SLradjb, sethfr30SLradjb};
var_set = {setnormadjb, setnormSLradjb};

%% Data from Sanghvi 
BSA = 1.8;
setSanghviD.tab = table();
setSanghviD.tab.HR = [78; 81];setSanghviD.tab.BPpp = [59;53];setSanghviD.tab.EDP = [8.2;16.2];
setSanghviD.tab.CO = [3.2;3.6]*BSA;
setSanghviD.tab.SV = [42;45]*BSA;
setSanghviD.tab.PCWP = setSanghviD.tab.EDP;
setSanghviD.tab.dPCWP = setSanghviD.tab.PCWP - setSanghviD.tab.PCWP(1);
setSanghviD.tab.SVrn = setSanghviD.tab.SV/setSanghviD.tab.SV(1);
setSanghviD.glyph = 'd';
setSanghviD.title = 'Sanghvi Dextran';

setSanghviI.tab = table();
setSanghviI.tab.HR = [85; 88];setSanghviI.tab.BPpp = [55.5;66.8];setSanghviI.tab.EDP = [8;18.4];
setSanghviI.tab.CO = [3.2;4]*BSA;
setSanghviI.tab.SV = [38;46]*BSA;
setSanghviI.tab.PCWP = setSanghviI.tab.EDP;
setSanghviI.tab.dPCWP = setSanghviI.tab.PCWP - setSanghviI.tab.PCWP(1);
setSanghviI.tab.SVrn = setSanghviD.tab.SV/setSanghviD.tab.SV(1);
setSanghviI.glyph = 'd';
setSanghviI.title = 'Sanghvi Infusion';

var_set = {setnormb, setnormSLradjb, setSanghviD, setSanghviI};
%% data from London and Safar 1989
BSA = 1.8;

setLondon.tab = readtable("..\..\Data\LondonSafar1989.csv");
setLondon.tab.SV = setLondon.tab.SI*BSA;
setLondon.tab.CO = setLondon.tab.CI*BSA;
setLondon.tab.dPCWP = setLondon.tab.PCWP - setLondon.tab.PCWP(1);
setLondon.tab.SVrn = setLondon.tab.SV/setLondon.tab.SV(1);

setLondon.glyph = '^';
setLondon.title = 'London1989';

var_set = {setnormb, setnormSLradjb, setSanghviD, setSanghviI, setLondon};
i_baseline = find(min(setnormSLradjb.tab.dVol.^2) == setnormSLradjb.tab.dVol.^2);
PCWP_target = setnormSLradjb.tab.PCWP(i_baseline) + setLondon.tab.dPCWP(end);
SV_target = setnormSLradjb.tab.SV(i_baseline)*setLondon.tab.SVrn(end);

disp("SetnormSLradjb BASE> SV=" + setnormSLradjb.tab.SV(i_baseline) + " @PCWP " + setnormSLradjb.tab.PCWP(i_baseline) + " mmHg");
disp("SetnormSLradjb target> SV=" + SV_target + " @PCWP " + PCWP_target + " mmHg");

%% data from Freitas1965
setFreitas.tab = readtable("..\..\Data\Freitas1965.csv");
setFreitas.tab.SV = setFreitas.tab.SI*BSA;
setFreitas.tab.CO = setFreitas.tab.CI*BSA;
setFreitas.tab.PCWP = setFreitas.tab.PLA;
setFreitas.tab.EDP = setFreitas.tab.PLA;
setFreitas.tab.dPCWP = setFreitas.tab.PCWP - setFreitas.tab.PCWP(2); % the first data point is venous pooling
setFreitas.tab.SVrn = setFreitas.tab.SV/setFreitas.tab.SV(2); % the first data point is venous pooling
setFreitas.glyph = '<';
setFreitas.title = 'Freitas1965';
var_set = {setnormb, setnormSLradj, setnormSLrOpt, setSanghviD, setSanghviI, setLondon, setFreitas};
checkAllVarSets(var_set)   

%% Data Bngaard-Nielsen2010
% Volume unloading by LBNP at normal
setBungaard.tab = table();
setBungaard.tab.PCWP = [2.398983482;4.899618806;9.758576874];
setBungaard.tab.SV = [73.65957447;93.31914894;109.7021277];
setBungaard.tab.dPCWP = setBungaard.tab.PCWP - setBungaard.tab.PCWP(3); % the first two data points are venous pooling
setBungaard.tab.SVrn = setBungaard.tab.SV/setBungaard.tab.SV(3); % the first two data point is venous pooling
setBungaard.glyph = '>';
setBungaard.title = 'Bungaard2010';
%%
var_set = {setnormb,setnormSLrOpt, setnormSLrOptb, setnormSLrOptnb, setOlsonsnb, setSanghviD, setSanghviI, setLondon, setFreitas};
checkAllVarSets(var_set)   


%%
close all;
cocomp = 5.7;
bpmcomp = 98;
xl = xlim;
yl = ylim;
% normal
plot([0, 0], [yl(1), yl(2)], 'k--');
%  plot([xl(1), xl(2)], [cocomp, cocomp], 'k--');
plot([xl(1), xl(2)], [bpmcomp, bpmcomp], 'k--');

plotHFVol_draw(var_set, 'dVol', 'CO', prenom);
plotHFVol_draw(var_set, 'dVol', 'EDV', prenom);
plotHFVol_draw(var_set, 'dVol', 'CVP', prenom);
plotHFVol_draw(var_set, 'dVol', 'PCWP', prenom);
plotHFVol_draw(var_set, 'dVol', 'SV', prenom);
plotHFVol_draw(var_set, 'dVol', 'EF', prenom);
plotHFVol_draw(var_set, 'dVol', 'BPm', prenom);
plotHFVol_draw(var_set, 'CO', 'CVP', prenom);
plotHFVol_draw(var_set, 'CO', 'PCWP', prenom);
plotHFVol_draw(var_set, 'CO', 'EDP', prenom);
plotHFVol_draw(var_set, 'CO', {'PCWP', 'EDP', 'ESP'}, prenom, {'*:', 'd-', '--'});
plotHFVol_draw(var_set, 'dPCWP', 'CO', prenom);
plotHFVol_draw(var_set, 'EDP', 'CO', prenom);
plotHFVol_draw(var_set, 'dPCWP', 'SVrn', prenom);
xlim([-10, 20])
ylim([0.2, 1.4])
ylim([40, 180])
xlim([0, 250])
plotHFVol_draw(var_set, 'dVol', 'EF', prenom);
plotHFVol_draw(var_set, 'dVol', 'BPm', prenom);
plotHFVol_draw(var_set, 'dVol', 'eGFR', prenom);
plotHFVol_draw(var_set, 'dVol', 'V_Isf_Exc', prenom);
plotHFVol_draw(var_set, 'dVol', 'p_int', prenom);
plotHFVol_draw(var_set, 'EDP', 'CO', prenom);
plotHFVol_draw(var_set, 'EDP', 'SV', prenom);
plotHFVol_draw(var_set, 'HR', 'EDP', prenom);

plotHFVol_draw(var_set, 'dPCWP', 'BPm', prenom);
plotHFVol_draw(var_set, 'dPCWP', 'HR', prenom);
plotHFVol_draw(var_set, 'dVol', 'HR', prenom);
plotHFVol_draw(var_set, 'CO', {'BPs', 'BPd'}, prenom);
plotHFVol_draw({setnormSLradj}, 'dVol', 'EDV', prenom);
plotHFVol_draw(var_set, 'HR', 'CO', prenom);
plotHFVol_draw(var_set, 'CVP', 'CO', prenom);
plotHFVol_draw(var_set, 'PCWP', 'CO', prenom);
plotHFVol_draw(var_set, 'CO', 'PCWP', prenom);
plot(3:10, (3:10)*1.2, 'k--');
plot(3:10, (3:10)*3.4, 'k--');
co = 0:7;pcwp = 1.2*co;plot(pcwp, co);

plotHFVol_draw(var_set, 'CO', 'BPs', prenom);

plotHFVol_draw(var_set, 'dVol', 'VP_r', prenom);

plotHFVol_draw(var_set, 'dVol', 'Pwr_r', prenom);
plotHFVol_draw(var_set, 'CO', {'Pwr_LV', 'Pwr_RV'}, prenom, {'+-' ,'x:'});
plotHFVol_draw(var_set, 'CO', {'Pwr_LVr', 'Pwr_RVr'}, prenom, {'+-' ,'x:'});
plotHFVol_draw(var_set, 'dVol', {'Pwr_LVr', 'Pwr_RVr'}, prenom, {'+-' ,'x:'});
plotHFVol_draw(var_set, 'PCWP', {'Pwr_LVr', 'Pwr_RVr'}, prenom, {'+-' ,'x:'});
plotHFVol_draw(var_set, 'CO', {'Pwr_LVrn', 'Pwr_RVrn'}, prenom, {'+-' ,'x:'});
plotHFVol_draw(var_set, 'dVol', 'Pwr_LV', 'on Exercise level');
plotHFVol_draw(var_set, 'dVol', 'Pwr_RV', 'on Exercise level');
plotHFVol_draw(var_set, 'CO', 'Pwr_RV', 'on Exercise level');
plotHFVol_draw(var_set, 'dVol', 'Pwr_LVr', prenom);
plotHFVol_draw(var_set, 'dVol', 'Pwr_RVr', prenom);


plotHFVol_draw(var_set, 'dVol', {'Pwr_LV', 'Pwr_RV'} , '');
plotHFVol_draw({setnormadj, sethfr25adj}, 'dVol', {'Pwr_LV', 'Pwr_RV'} , ' normal with baroreflex');
%% Comparison with LVSWI data
plotHFVol_draw(var_set, 'EDP', 'LVSWI', prenom);
plot(8, 46, 's'); % sanghvi
plot(18, 52, 's'); % sanghvi
plot(8, 49, 'd'); % sanghvi
plot(16, 52, 'd'); % sanghvi

%% ventricular energetcs siomulated
plotHFVol_draw(var_set, 'vol', 'PCWP', prenom);
plotHFVol_draw(var_set, 'vol', 'HR', prenom);
plotHFVol_draw(var_set, 'PCWP', {'Pwr_LVr', 'Pwr_RVr'}, '(Ratio of normal)');
plotHFVol_draw(var_set, 'PCWP', {'Pwr_LV', 'Pwr_RV'}, '(W)');

%% ventricular energetics
% Compensated by CO at 0.25 and 2.4
pwrlv = [1.25, 1.35, 1.2];
pwrrv = [0.25, 0.3, 0.48];
pwrlv2 = [1.25, 1.32, 1.12];
pwrrv2 = [0.25, 0.28, 0.39];

figure('Name', 'Ventricular power');hold on;
plot(1:size(pwrlv, 2), pwrlv, 'b*-');
plot(1:size(pwrrv, 2), pwrrv, 'r*-');
plot(1:size(pwrlv2, 2), pwrlv2, 'b+:');
plot(1:size(pwrrv2, 2), pwrrv2, 'r+:');
title('Ventricular power output');
set(gca,'XTick',1:3,'XTickLabel',{'Norm, EF 60%', 'EF 50%, +250ml', 'EF 40%, +2.4L'});
ylabel('Power output (W)');
legend('LV power (CO compensated)', 'RV power (CO compensated)','LV power (BPm compensated)','RV power (BPm compensated)','Location', 'Best');

figure('Name', 'Ventricular power (% of baseline)');hold on;
plot(1:size(pwrlv, 2), pwrlv/pwrlv(1)*100, 'b*-');
plot(1:size(pwrrv, 2), pwrrv/pwrrv(1)*100, 'r*-');
plot(1:size(pwrlv2, 2), pwrlv2/pwrlv2(1)*100, 'b+:');
plot(1:size(pwrrv2, 2), pwrrv2/pwrrv2(1)*100, 'r+:');
title('Ventricular power output (% of baseline)');
set(gca,'XTick',1:3,'XTickLabel',{'Norm, EF 60%', 'EF 50%, +250ml', 'EF 40%, +2.4L'});
ylabel('Power output (W)');
legend('LV power (CO compensated)', 'RV power (CO compensated)','LV power (BPm compensated)','RV power (BPm compensated)','Location', 'Best');
%% Simplify the Baroreflex 
var_set = {setnormadjb, setnormSLradjb}
pps = [];
hrs = [];
phis = []
cars = []
aors = []
for i = 1:size(var_set,2)
%     var_set{i}.glyph = 'kx';
    pps = [pps; var_set{i}.tab.BPpp];
    hrs = [hrs; var_set{i}.tab.HR];
    cars = [cars; var_set{i}.tab.dCar_d];
    aors = [aors; var_set{i}.tab.dAor_d];    
end
phis = (hrs/60 / 1.06667 - 1 )/2.625 + 0.25;
    
plotHFVol_draw(var_set, 'BPpp', 'HR', prenom);


ppss = 5:1:70;
phiEst = 1.518 * ppss.^(-0.5048);
hrsEst = 1.06667*(1 + 2.625*(phiEst - 0.25));

figure;clf;hold on;
% plot(ppss, phiEst, '-', 'LineWidth', 1);
plot(pps, phis, 'x');
plot([ppss(1), ppss(end)], [0.25, 0.25], 'k--');
xlabel('brachial pulse pressure');
ylabel('\Phi');
legend('Simulated', 'Simplified', 'Baseline');

figure;clf;hold on;
% plot(ppss, hrsEst*60, '-', 'LineWidth', 1);
plot(pps, hrs, '+', 'MarkerSize', 12);
plot([ppss(1), ppss(end)], [64, 64], 'k--');
xlabel('brachial pulse pressure');
ylabel('HR');
legend('Simplified', 'Simulated', 'Baseline');
title('Long-term baroreflex resetting to steady state')

ddist = 0.01: 0.01:0.3;
phiEst = (0.04637) ./ (ddist + 0.05604);
hrsEst = 1.06667*(1 + 2.625*(phiEst - 0.25));
       
figure;
clf;hold on;
plot(ddist, hrsEst*60, '-', 'LineWidth', 1);
plot(cars, hrs, 'b+', 'MarkerSize', 12);
plot(aors, hrs, 'r+', 'MarkerSize', 12);
% plot([ppss(1), ppss(end)], [64, 64], 'k--');
xlabel('Vessel distension difference');
ylabel('HR');
% legend('Simplified', 'Simulated', 'Baseline');
title('Long-term baroreflex resetting to steady state');


% set_norm = importVolLoadResultFile("..\..\Results\nobaro_VolLoad_result_normal_.csv");
% set_hfr = importVolLoadResultFile("..\..\Results\nobaro_VolLoad_result_hfref_.csv");
% set_hfp = importVolLoadResultFile("..\..\Results\nobaro_VolLoad_result_hfpef_.csv");
% 
% prenom = 'NO baroreflex';
% var_set = {set_norm, set_hfr, set_hfp};
% 
% plotHFVol_draw(var_set, 'dVol', 'eGFR', prenom);
% plotHFVol_draw(var_set, 'dVol', 'CO', prenom);
% plotHFVol_draw(var_set, 'dVol', 'HR', prenom);
% plotHFVol_draw(var_set, 'HR', 'CO', prenom);
% plotHFVol_draw(var_set, 'CVP', 'CO', prenom);
% plotHFVol_draw(var_set, 'PCWP', 'CO', prenom);
% plotHFVol_draw(var_set, 'CVP', 'PCWP', prenom);
% plotHFVol_draw(var_set, 'CO', 'BPs', prenom);