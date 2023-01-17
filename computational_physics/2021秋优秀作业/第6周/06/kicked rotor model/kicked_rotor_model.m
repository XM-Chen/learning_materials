clear 
clc

% 实现 Kicked totor model的模拟
% 迭代法
% 只绘制一种情况
% 改变x0和k的程序在kicked_rotor_model_k.m与kicked_rotor_model_x0.m中

% 预设
k=1.5;


% 初始条件
x0=2.6;
p0=1;

% 迭代次数
N=1e5;

% 初始化
x=zeros(1,N+1);
p=zeros(1,N+1);

% 初始条件
% 对2pi取模，并且保留正负号
x(1)=mod(abs(x0),2*pi)*sign(x0);
p(1)=mod(abs(p0),2*pi)*sign(p0);



% 迭代
for i=1:N
    % p迭代
    p(i+1)=p(i)+k*sin(x(i));
    
    % 模2*pi
    p(i+1)=mod(abs(p(i+1)),2*pi)*sign(p(i+1));
    
    % x迭代
    x(i+1)=x(i)+p(i+1);
    
    % 模2*pi
    x(i+1)=mod(abs(x(i+1)),2*pi)*sign(x(i+1));
end



% 绘制相空间图像
figure

% 全屏幕显示
set(gcf,'outerposition',get(0,'screensize'));

% 绘制散点图
scatter(x,p,'.b');

% 坐标轴设置
xlabel('x');
ylabel('p');
legend('k=1.5');
title('kicked rotor model')

