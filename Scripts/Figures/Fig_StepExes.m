% plot figure 11
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
path = '../../Results/';
filenames = ["ExStep_base", "ExStep_chronotropy", "ExStep_inotropy", "ExStep_noVc", "ExStep_noAR", "ExStep_Tilted60", "ExStep_Tilted60_noVC", "ExStep_Tilted60_noVCLin"];
tities = ["A: Supine normal", "B: No chronotropy", "C: No inotropy", "D: No VC", "E: No arteriolar vasodilation", "F: 60° HUT", "G: 60° HUT, no VC", "H: 60° HUT, no VC, lin PV"];
color_b = [28, 108, 200]/255;
color_r = [238, 46, 47]/255;
color_g = [0, 140, 72]/255;
color_m = [226, 113, 199]/255;
color_lb = [182, 226, 255]/255;
%% covenrt to sdfs
% for i = 2:size(filenames, 2)
%     filename = filenames(i);
%     disp('Filename ' + filename + '...')
%     system('"c:\Program Files\Dymola 2021x\bin\dsres2sdf.exe" ../../Results/' + filename + '.mat ' + filename + '.sdf');
%     pbm_avst = h5read(filename + '.sdf', '/brachial_pressure_mean')/mmHg2SI;
%    
%     figure();title(filename);plot(pbm_avst)
% end
% %%
% pbm_avst = h5read(filenames(1) + '.sdf', '/brachial_pressure_mean')/mmHg2SI;
% plot(pbm_avst)
%%
fig = figure(11);clf;

tw = 17;
th = 12;
set(fig, 'Units', 'Centimeters', 'Position', [2, 4, tw, th])

h_sp = 0.01
v_sp = 0.06

h_ins = 0.06
v_ins = 0.06

h_s = (1-2*h_ins - h_sp*3)/3
v_s = (1 - 1*v_ins - v_sp*3)/3


i = 1;
pow_lv = {}
for row = 3:-1:1
     for col = 1:3
         
        hpos = h_ins + (col-1)*(h_s + h_sp);
        vpos = v_ins + (v_sp/2 + (row - 1)*(h_s + h_sp));
        ax = axes
        ax.Position = [hpos, vpos, h_s, v_s];
        
         if i > size(filenames)
             %% final figure
             hpos = h_ins + (col-1)*(h_s + h_sp);
             vpos = v_ins + (v_sp/2 + (row - 1)*(h_s + h_sp));

             hold on;
             title('I: LV power output')
             lege = {};
             plotEner = [1 2 3 6 7];
             colors = [color_b; color_r; color_g; color_m; color_lb]
             yyaxis left;
             set(gca,'ytick',[]);
             yyaxis right;
             for j = 1:size(plotEner, 2)
                 plot([0:10:100], pow_lv{plotEner(j)}, 'LineWidth', 1, 'Color', colors(j, :));
                 lege{plotEner(j)} = char(tities{plotEner(j)}(1));
             end
             g = gca; g.YAxis(2).Color = [0, 0, 0];
             xlabel('Exercise (% of max)');
             ylabel('Power output (W)');
             
            leg = legend(lege{plotEner}, 'Location', 'NorthWest', 'Orientation', 'Vertical');
            leg.ItemTokenSize = [10, 150];
            ax.Position = [hpos + 0.05, vpos, h_s - 0.05, v_s];

            %%
            continue
         end
         
        
        s = plotExSubplot(path + filenames(i) + '.mat', tities(i), false, i == 3);
        if col == 1 
            % disable yy axis
            yyaxis right;
            ylabel('');
            set(gca,'ytickLabel',[]);
        elseif col == 2 & row ~= 1
            % disable both axis
            yyaxis left;
            ylabel('');
            set(gca,'ytickLabel',[]);
            yyaxis right;
            ylabel('');
            set(gca,'ytickLabel',[]);
        elseif col == 3 | (col == 2 & row == 1)
            % disable y axis
            yyaxis left;
            ylabel('');
            set(gca,'ytickLabel',[]);
        end
        
        if row ~= 1
            xlabel('');
            set(gca,'xtickLabel',[]);
        end
        
            
          
            
        pow_lv{i} = s.pow_lv;
        i = i + 1;
     end
end
%%
exportgraphics(gcf,'fig_R_StepExImp.png','Resolution',200)
exportgraphics(gcf,'fig_R_StepExImp.pdf', 'ContentType','vector')

