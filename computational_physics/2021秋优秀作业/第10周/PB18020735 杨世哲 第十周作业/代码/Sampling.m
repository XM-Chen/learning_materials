%Sampling
format long;
%% Sampling 1
clear;
N = 1e6;
NormFac = 0.7.*sqrt(pi./10) + 0.9.*sqrt(3/8).*pi.*exp(9/8).*(besseli(-0.25, 1.125) + besseli(0.25, 1.125));
p1 = @(x) 0.7.*exp(-10.*(x+10).^2) + 0.9.*exp(-x.^4 + 3.*x.^2);
Dist1 = @(x) (0.7.*exp(-10.*(x+10).^2) + 0.9.*exp(-x.^4 + 3.*x.^2))./NormFac;
New1 = @(~) 4.5.*(rand - 0.5) - 5 + 5.*sign(rand - 0.5);
xi1 = StableMetropolis(N, 1, p1, New1, Dist1);
%Calculate distribution
[x1, fx1] = Distribution(xi1, 200);
fig1 = figure;
scatter(x1{:, :}, fx1, 100, '.', 'DisplayName', 'Data');
hold on;
y1 = p1(x1{:, :});
plot(x1{:, :}, y1, 'LineWidth', 1.5, 'DisplayName', 'Aim Distribution');
legend;
savefig(fig1, 'Expotional.fig');
saveas(fig1, 'Expotional.png');

%% Sampling 2
clear;
N = 1e6;
m = 1;
a = 5;
b = 10;
epsilon = 0.1;
p2 = @(x) 1./(x(:, 1).^2 + m.^2) + 1./(x(:, 2).^2);
New2 = @(~) [5.*rand, 9.9.*rand + 0.1];
Fitting1 = @(~, ~) 1;
A = a./epsilon - a./b + atan(a./m);
Dist2 = @(x, y) (1./(x.^2 + m.^2) + 1./(y.^2))./A;
xi2 = StableMetropolis(N, 2, p2, New2, Dist2);
%Calculate distribution
[x2, fx2] = Distribution(xi2, [100, 100]);
y2 = Dist2(x2{:, :});
RelErr21 = (fx2 - y2)./y2;
fig21 = figure;
imagesc(x2{:}, RelErr21);
colorbar;
xlabel('x');
ylabel('y');
savefig(fig21, 'UnidentifiedDist1.fig');
saveas(fig21, 'UnidentifiedDist1.png');
RelErr22 = (fx2 - y2)./max(y2, [], 'all');
fig22 = figure;
imagesc(x2{:}, RelErr22);
colorbar;
xlabel('x');
ylabel('y');
savefig(fig22, 'UnidentifiedDist2.fig');
saveas(fig22, 'UnidentifiedDist2.png');




