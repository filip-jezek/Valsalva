% read and plot valsalva data from Kosinksi et al. 2018
data1 = xlsread('..\..\data\Valsalva\norm004_BP.xlsx', 'A3:C98711');
data2 = xlsread('..\..\data\Valsalva\norm004_HR.xlsx', 'A5:M580');
t3 = xlsread('..\..\data\Valsalva\norm004_Pexp.xlsx', 'B2:JN2');
tps = xlsread('..\..\data\Valsalva\norm004_Pexp.xlsx', 'B6:JN15');
%%
t1 = data1(:, 1);
finap = data1(:, 2);
reBap = data1(:, 3);

t2 = data2(:, 1);
hr = data2(:, 8);
lvet = data2(:, 10);
sv = data2(:, 11);
tpr = data2(:, 12);

tp = mean(tps, 1);

%%
figure(1);clf;hold on;
% plot(t1, finap);
plot(t1, reBap);

plot(t2, hr);
plot(t2, lvet);
plot(t2, sv);
plot(t2, tpr);
legend('BP', 'HR', 'LVET', 'SV', 'tpr')
%%
figure(1);clf;
t3_ = (1:size(tp, 2))*1.0 - 67.0;
plot(t3_, tp, 'r*-')


