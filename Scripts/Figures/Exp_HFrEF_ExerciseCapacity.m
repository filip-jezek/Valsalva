% plot figure 11
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
path = '../../Results/';
filenames = ["ExStep_base", "CVS_ExCap_Corvia","CVS_ExCap_Corvia_IncrRPulm","CVS_ExCap_Corvia_DecRPulm", "CVS_HFrEF_15_ExerciseCapacity", "CVS_HFrEF_15_ExerciseCapacity_Scaled","CVS_HFrEF_15_ExCap_Sc_Corvia","CVS_HFrEF_15_ExCap_Sc_Corvia_IncRPulm","CVS_HFrEF_15_ExCap_Sc_Corvia_DecRPulm"];
titIes = filenames;
% tities = ["A: HFrEF 15", "B: HFrEF 15 scaled to start at 90bpm", "C: No inotropy", "D: No VC", "E: No arteriolar vasodilation", "F: 60° HUT", "G: 60° HUT, no VC", "H: 60° HUT, no VC, lin PV"];
color_schema;


for i = 1:size(filenames, 2)
    figure(i);clf;hold on;
    s = plotExSubplot(path + filenames(i) + '.mat', titIes(i), false, false);
    yyaxis left;ylim([0, 250]);
    plot(s.t_ax, s.pcwp, '.--', 'Color', color_s);    
    plot(s.t_ax, s.ef.*100, '*-', 'Color', color_r);    
    yyaxis right;ylim([0, 25]);
    plot(s.t_ax, s.q_c, '+-','Color',color_g);
end;