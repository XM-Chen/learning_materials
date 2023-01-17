%Integration
format long;
%% Integral 1
clear;
%--------------Parameters and distributions--------------%
N = 1e5;
NormFac = pi.^2./24;
p = @(x) exp(- x(:, 1).^2 - 4.*x(:, 2).^2 - 9.*x(:, 3).^2 - 16.*x(:, 4).^2);
Dist = @(x1, x2, x3, x4) exp(- x1.^2 - 4.*x2.^2 - 9.*x3.^2 - 16.*x4.^2)./NormFac;
New = @(~) [7.*(rand - 0.5), 3.5.*(rand - 0.5), 2.4.*(rand - 0.5), 1.75.*(rand - 0.5)];
%--------------Sampling--------------%
x = StableMetropolis(N, 4, p, New, Dist);
y = x(:, 1).^2.*exp(x(:, 1).*x(:, 2) + x(:, 1).*x(:, 3) + x(:, 1).*x(:, 4) + x(:, 2).*x(:, 3) + x(:, 2).*x(:, 4) + x(:, 3).*x(:, 4));
%--------------Calculate integral and plot figures of distribution of x1~x4--------------%
fig = figure;
hold on;
for k = 1:4
    [xk, fxk] = Distribution(x(:, k), 200);
    A = GaussianFitting(xk{:}, fxk);
    fprintf('sqrt(2).*%d.*Sigma_%d = %f\n', k, k, A(1).*sqrt(2).*k);
    lgd = strcat('x_', num2str(k));
    scatter(xk{:}, fxk, 100, '.', 'DisplayName', lgd);
end
legend;
savefig(fig, 'x1-x4Distribution.fig');
saveas(fig, 'x1-x4Distribution.png');
Iz = mean(y)./96;
fprintf('I_0 = %e\n', Iz);

%This integral could be calculated analytically, so there is no need to
%iterate to check the accuracy.

% % % %--------------Iterate to make the result stable--------------%
% % % MacIter = 10;
% % % IterTime = 1;
% % % ErrMax = 1e-2;
% % % N0 = numel(x0);
% % % x1 = Metropolis(N0, p, New)./sqrt(2);
% % % x2 = Metropolis(N0, p, New)./2./sqrt(2);
% % % x3 = Metropolis(N0, p, New)./3./sqrt(2);
% % % x4 = Metropolis(N0, p, New)./4./sqrt(2);
% % % y = x1.^2.*exp(x1.*x2 + x1.*x3 + x1.*x4 + x2.*x3 + x2.*x4 + x3.*x4);
% % % Iz_incr = mean(y).*pi.^2./24;
% % % Iz = (Iz.*IterTime + Iz_incr)./(IterTime + 1);
% % % fprintf('I_%d = %f\n', IterTime, Iz);
% % % err = abs(1 - Iz_incr./Iz);
% % % while err > ErrMax
% % %     IterTime = IterTime + 1;
% % %     if IterTime > MaxIter
% % %         warning('Unable to calculate accurate integral.\n');
% % %         break;
% % %     end
% % %     x1 = Metropolis(N0, p, New)./sqrt(2);
% % %     x2 = Metropolis(N0, p, New)./2./sqrt(2);
% % %     x3 = Metropolis(N0, p, New)./3./sqrt(2);
% % %     x4 = Metropolis(N0, p, New)./4./sqrt(2);
% % %     y = x1.^2.*exp(x1.*x2 + x1.*x3 + x1.*x4 + x2.*x3 + x2.*x4 + x3.*x4);
% % %     Iz_incr = mean(y).*pi.^2./24;
% % %     Iz = (Iz.*IterTime + Iz_incr)./(IterTime + 1);
% % %     fprintf('I_%d = %f\n', IterTime, Iz);
% % %     err = abs(1 - Iz_incr./Iz);
% % % end
% % % fprintf('I = %f\n', Iz);
% % % fprintf('RelError ~ %f\n', err);

%% Integral 2
clear;
%--------------Parameters--------------%
N = 1e6;
Lambda = 1;
s = 1.5;
Mu = 0.1;
% p = 0.8;
P = [p 0 0];
%--------------Distribution function--------------%
pk = @(k) k.^2./(k.^2 + Mu.^2);
Newk = @(~) Lambda./s + rand.*Lambda.*(1 - 1./s);
Ak = Lambda - Mu.*atan(Lambda./Mu) - Lambda./s + Mu.*atan(Lambda./Mu./s);
Distk = @(k) k.^2./(k.^2 + Mu.^2)./Ak;
Theta = @(theta) sin(theta);
NewTheta = @(~) rand.*pi;
DistTheta = @(theta) 0.5.*sin(theta);
%--------------Sample points--------------%
K = StableMetropolis(N, 1, pk, Newk, Distk);
theta = StableMetropolis(N, 1, Theta, NewTheta, DistTheta);
N1 = numel(K);
N2 = numel(theta);
if N1 > N2
    theta = [theta; Metropolis(N1 - N2, 1, Theta, NewTheta, theta(end))];
elseif N1 < N2
    K = [K; Metropolis(N2 - N1, 1, pk, Newk, K(end))];
end
N0 = max(N1, N2);
phi = 2.*pi.*rand(size(K));%phi
I = IntegrateKT(K, theta, phi, N0, Lambda, Mu, s, P);
%--------------Iterate to make the result stable--------------%
MacIter = 10;
IterTime = 1;
ErrMax = 1e-2;
K = Metropolis(N0, 1, pk, Newk, K(end));
theta = Metropolis(N0, 1, Theta, NewTheta, theta(end));
phi = 2.*pi.*rand(size(K));
I_incr = IntegrateKT(K, theta, phi, N0, Lambda, Mu, s, P);
I1 = (I.*IterTime + I_incr)./(IterTime + 1);
% fprintf('I_%d = %f\n', IterTime, I1);
err = abs(1 - I./I1);
I = I1;
while err > ErrMax
    IterTime = IterTime + 1;
    if IterTime > MaxIter
        warning('Unable to calculate accurate integral.\n');
        break;
    end
    K = Metropolis(N0, 1, pk, Newk, K(end));
    theta = Metropolis(N0, 1, Theta, NewTheta, theta(end));
    phi = 2.*pi.*rand(size(K));
    I_incr = IntegrateKT(K, theta, phi, N0, Lambda, Mu, s, P);
    I1 = (I.*IterTime + I_incr)./(IterTime + 1);
%     fprintf('I_%d = %f\n', IterTime, I1);
    err = abs(1 - I./I1);
    I = I1;
end
fprintf('IterTime = %d     p = %e     I = %e     RelError ~ %e\n', IterTime, p, I, err);
% fprintf('RelError ~ %f\n', err);

%--------------Function for calculating the integral--------------%
function I = IntegrateKT(K, theta, phi, N, Lambda, Mu, s, P)
    %Divide
    k1 = K(1:N/2);
    theta1 = theta(1:N/2);
    phi1 = phi(1:N/2);
    k1x = k1.*sin(theta1).*cos(phi1);
    k1y = k1.*sin(theta1).*sin(phi1);
    k1z = k1.*cos(theta1);
    k2 = K(N/2+1:end);
    theta2 = theta(N/2+1:end);
    phi2 = phi(N/2+1:end);
    k2x = k2.*sin(theta2).*cos(phi2);
    k2y = k2.*sin(theta2).*sin(phi2);
    k2z = k2.*cos(theta2);
    %Integral
    I = (Lambda - Mu.*atan(Lambda./Mu) - Lambda./s + Mu.*atan(Lambda./Mu./s)).^2./(4.*pi.^4)./mean(sum((P - [k1x k1y k1z] - [k2x k2y k2z]).^2 + Mu.^2, 2));
end

% % % % %--------------Function for fitting the distribution of k--------------%
% % % % function para = Kfitting(x, y)
% % % %     %Reshape x & y to column vectors.
% % % %     x = reshape(x, numel(x), 1);
% % % %     y = reshape(y, numel(y), 1);
% % % %     z = 1./y;
% % % %     %Fitting
% % % %     f = @(x) [x.^2 ones(size(x))];
% % % %     para = OLSfitting(f, x, z);
% % % % end

% % % % %--------------Function for fitting the distribution of theta--------------%
% % % % function para = ThetaFitting(x, y)
% % % %     %Reshape x & y to column vectors.
% % % %     x = reshape(x, numel(x), 1);
% % % %     y = reshape(y, numel(y), 1);
% % % %     %Fitting
% % % %     f = @(x) [sin(x) cos(x) ones(size(x))];
% % % %     para = OLSfitting(f, x, y);
% % % % end
