tic
clear;clc;
v0 = 0;
t_end = 10;
t_step = 1000;
dt = t_end / t_step;
v = zeros(t_step + 1, 1);
for i = 2:t_step + 1
    v(i) = v(i - 1) + (-(1 + (i - 1) * dt) * v(i - 1) + (4 * rand - 2)) * dt;
end