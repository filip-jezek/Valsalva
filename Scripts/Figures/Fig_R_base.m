%% plots the baseline results
datafile = '../../Results2/CardiovascularSystem.mat'
% import the dymload util
addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')
dl = dymload(datafile)
%%
mmHg2SI = 133.322;
ml2SI = 1e-6;
%%
time = dymget(dl, 'Time');
t_interval = [58.4, 60]; % interval in seconds
td = t_interval(2) - t_interval(1)
% i_int = [find(t >= t_interval(1), 1), find(t >= t_interval(2), 1)-1];
i_int = (time >= t_interval(1) & time <= t_interval(2));

t = time - t_interval(1);
pb = dymget(dl, 'brachial_pressure')/mmHg2SI;
psa = dymget(dl, 'heartComponent.sa.pressure')/mmHg2SI;
plv = dymget(dl, 'heartComponent.mitralValve.q_out.pressure')/mmHg2SI;
pla = dymget(dl, 'heartComponent.mitralValve.q_in.pressure')/mmHg2SI;

vlv = dymget(dl, 'V_LV')/ml2SI;
vla = dymget(dl, 'V_la')/ml2SI;

ppa = dymget(dl, 'heartComponent.pa.pressure')/mmHg2SI;
prv = dymget(dl, 'heartComponent.tricuspidValve.q_out.pressure')/mmHg2SI;
pra = dymget(dl, 'heartComponent.tricuspidValve.q_in.pressure')/mmHg2SI;

vrv = dymget(dl, 'heartComponent.ventricles.V_RV')/ml2SI;
vra = dymget(dl, 'heartComponent.ra.volume')/ml2SI;

qmv = dymget(dl, 'q_mv')/ml2SI/1000*60;
%%
fig1 = figure(1);clf;
set(gcf, 'DefaultAxesFontSize', 10, 'defaultLineLineWidth',1.0);

s_a1 = subplot(4, 2, 1);
% a_a1 = axes()
% a_a1 = gca()

hold on;
% title('A: Left ventricular pressures and volumes', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
title('A: Left ventricular pressures and volumes');
plot(t(i_int), plv(i_int), 'b', 'LineWidth', 1);
plot(t(i_int), psa(i_int), 'r', 'LineWidth', 1);
plot(t(i_int), pla(i_int), 'm', 'LineWidth', 1);
plot(t(i_int), pb(i_int));
set(gca,'xtick',[])
leg = legend('P LV', 'P Asc Aor', 'P brach art', 'P LA', 'Location', 'SouthWest')
leg.ItemTokenSize = [10, 2];
% legend('boxoff')
ylim([0, 120]);
ylabel('Pressure (mmHg)')
s_a1.Clipping = 'off';
xlim([0 td])
%%    
% volumes
s_a2 = subplot(4, 2, 3);hold on;
plot(t(i_int), vlv(i_int), 'b');
plot(t(i_int), vla(i_int), 'r');
set(gca,'xtickMode', 'auto')
ylim([0, 150]);
xlim([0 td])
xlabel('t (s)');
ylabel('Volume (ml)')
% pos = get(s1, 'Position');
% set(s1, 'Position', [0.05, 0.5, 0.9,0.4])
% set(s2, 'Position', [0.05, 0.05, 0.9,0.4])

% B
s_b1 = subplot(4, 2, 2);cla;
hold on;
title('B: Right ventricular pressures and volumes');
plot(t(i_int), prv(i_int), 'b', 'LineWidth', 1);
plot(t(i_int), ppa(i_int), 'r', 'LineWidth', 1);
plot(t(i_int), pra(i_int), 'm', 'LineWidth', 1);
% set(gca,'xtick',[], 'ytick', [])
set(gca,'xtick',[])
leg = legend('P RV', 'P Pulm Ar', 'P RA')
leg.ItemTokenSize = [10, 2];
% ylim([0, 120])
xlim([0 td])
% volumes
s_b2 = subplot(4, 2, 4);hold on;
plot(t(i_int), vrv(i_int), 'b', 'LineWidth', 1);
plot(t(i_int), vra(i_int), 'r', 'LineWidth', 1);
set(gca,'xtickMode', 'auto', 'ytick', [])
ylim([0, 150]);
xlim([0 td])
leg = legend('V LV', 'V LA');
leg.ItemTokenSize = [10, 2];
xlabel('t (s)');
% pos = get(s1, 'Position');
% set(s1, 'Position', [0.05, 0.5, 0.9,0.4])
% set(s2, 'Position', [0.05, 0.05, 0.9,0.4])

% C
s_c = subplot(2, 2, 3);cla;hold on;
%%
% clf;hold on;
title('C: Mitral flow profile')


fill_pivot = find(i_int > 0, 1) + 65;
i_int_f1 = i_int(1:fill_pivot);
i_int_f2 = i_int;
i_int_f2(1:fill_pivot) = 0;
stop_pivot = fill_pivot + find(qmv(i_int_f2(fill_pivot+1:end)) < 0, 1);
i_int_f2(stop_pivot:end) = 0;

% fill([t(i_int_f1),fliplr(t(i_int_f1)),[t(fill_pivot); 0]], [qmv(i_int_f1), fliplr(qmv(i_int_f1));[0; 0]], 'b');
fill([t(i_int_f1); t(fill_pivot)], [qmv(i_int_f1); 0], [0.6 0.9 1]);
fill([t(fill_pivot);t(i_int_f2)], [0;qmv(i_int_f2)], [1 0.8 0.8]);

plot(t(i_int), qmv(i_int), 'Color', [0 0.4 1], 'LineWidth', 1);
% plot(t(i_int_f2), qmv(i_int_f2));
% fill([t(i_int), ])
xlabel('t (s)');
ylabel('Flow [L/min]');

%%
% D
s_d = subplot(2, 2, 4);cla;hold on;
title('D: Ventricular PV loop')
% plot(vra(i_int), pra(i_int));
% plot(vla(i_int), pla(i_int));
plot(vrv(i_int), prv(i_int), 'r', 'LineWidth', 1);
plot(vlv(i_int), plv(i_int), 'b', 'LineWidth', 1);
xlabel('Volume (ml)');
ylabel('Pressure (mmHg)');

%% align the axis
drawnow()
pos = get(s_a2, 'Position');

hm = 0.1;
vm = 0.05;
set(s_a1, 'Position', [hm, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
set(s_a2, 'Position', [hm, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])

set(s_b1, 'Position', [0.5+hm/2, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
set(s_b2, 'Position', [0.5+hm/2, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])

set(s_c, 'Position', [hm, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])
set(s_d, 'Position', [0.5 + hm/2, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])
% %%
% tw = 17;
% set(s_a1, 'Units', 'Centimeters', 'Position', tw*[hm, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
% set(s_a2, 'Units', 'Centimeters', 'Position', tw*[hm, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])
% 
% set(s_b1, 'Units', 'Centimeters', 'Position', tw*[0.5+hm/2, 0.75 + vm/2, (1 - 3*hm)/2,0.25 - vm])
% set(s_b2, 'Units', 'Centimeters', 'Position', tw*[0.5+hm/2, 0.5+vm, (1 - 3*hm)/2,0.25 - vm])
% 
% set(s_c, 'Units', 'Centimeters', 'Position', tw*[hm, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])
% set(s_d, 'Units', 'Centimeters', 'Position', tw*[0.5 + hm/2, 1.5*vm, (1 - 3*hm)/2,0.5 - 3*vm])
% 
% %% Draw boxes
% annotation(fig1,'rectangle', [0 0.5 0.5 0.5]);
% annotation(fig1,'rectangle', [0.5 0.5 0.5 0.5]);


%% save figure
tw = 17;
th = 17;
set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, tw, th], 'PaperUnits', 'Centimeters', 'PaperSize', [tw, th])
%%
exportgraphics(gcf,'fig_R_base.png','Resolution',300)
exportgraphics(gcf,'fig_R_base.pdf', 'ContentType','vector')
% print(gcf,'fig_R_base_300.png', '-dpng', '-r300')
% saveas(gcf,'fig1.svg')