function plotHFVol_draw(set, var1, var2, prenom, glyphs)
lw = 1;
if iscell(var2)
    var2_name = strjoin(var2(:), ' and ');
else
    var2_name = var2;
end

tit = [var1, ' vs. ', var2_name, ' ', prenom];

figure('Name', tit); hold on;
title(tit, 'Interpreter', 'None');
xlabel(var1, 'Interpreter', 'None');
ylabel(var2, 'Interpreter', 'None');
leg = {};
for i = 1:size(set, 2),
    s = set(i);
    if ~iscell(var2)
        plot(s{1}.tab.(var1), s{1}.tab.(var2), s{1}.glyph, 'LineWidth', lw);
        leg{end + 1} = s{1}.title;
    else
        for j = 1:size(var2, 2),
            leg{end + 1} = [s{1}.title , ' ', var2{j}];
            
            if nargin >= 5
                gt = glyphs(j);
                g = gt{1};
            elseif ~iscell(s{1}.glyph)
                g = s{1}.glyph;
%             else
%                 g = s{1}.glyph{j};
            end
            plot(s{1}.tab.(var1), s{1}.tab.(var2{j}), g, 'LineWidth', lw);
        end
    end            
end

legend(leg, 'Location', 'Best');