ReadFig1Data()
mmHg2SI = 133.322; % convert to Pa
dV2r = 1/100; % convert to 

p = GuytonCR65fig1.p*mmHg2SI; % pressure in Pa
Vr = GuytonCR65fig1.Vol*dV + 1.0 % Relative change in volume

p1 =   6.443e-09
p2 =   4.159e-05
p3 =       1.045
a =        1.06  
b =  -5.451e-05  
c =     -0.7703  
d =  -0.0007743  

t = (-3000:10:6000);
% x = p;
x = t;
y = max(0, p1*x.^2 + p2*x + p3) + max(0,a*exp(b*x) + c*exp(d*x));

figure(1);clf;hold on;
plot([min(p), max(p)], [1, 1], '--k');
plot(p/mmHg2SI , Vr, 'r*');
plot(x/mmHg2SI, y, '-b');
xlabel('Pressure [Pa]');
ylabel('Relative change in total weight [1]');
