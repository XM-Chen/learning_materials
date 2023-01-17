%Logistic map数值求解

clear
format long

%参数γ
r = 0.892;

%x的初值
x(1) = 0.3;
N = 50000;

%用迭代式进行迭代
for i = 1:N
x(i+1) = 4*r*x(i)*(1-x(i));
end

%作出迭代结果图
figure()
plot(x,".");
title("r ="+num2str(r))
pause(0.01)
ylabel('xn');
xlabel('n');

(r-0.25)/r
