tic
clear;clc;
N = 10;   % 循环次数
s = zeros(1, N);
for i = 1:N
    %% parameters
    tspan = [0 10];
    v0 = 0;
    
    %% Runge-Kutta
    [t, v] = ode45(@(t, v) odefun1(t, v), tspan, v0);
    s(i) = v(end);
    if i < N
        clear v;
        clear t;
    end
end
v_ave = sum(s) / N;
sigma2 = sum(s .^ 2) / N - v_ave ^ 2;
toc