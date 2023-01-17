clear 
clc

% 实现 Kicked totor model的模拟
% 迭代法
% 大概有点跑题
% 实现了粒子轨迹的追踪
% 具体在PDF中详细解释

% 预设
m=1;
k=5;
t0=1;

% 时间步长
dt=0.001;

% 迭代次数
N=1e5;

% 初态
x0=0;
p0=-8;

% 初始化，减少运算时间
x=zeros(1,N+1);
p=zeros(1,N+1);

% 初始条件
x(1)=x0;
p(1)=p0;

% 初始时刻接近0但不是0，便于判断是否到达动量突变点。
% 时间格点
t=dt/100:dt:dt/100+N*dt;


for i=1:N
    x(i+1)=x(i)+dt*p(i);
        
    % 判断是否到达动量突变点
    if abs(mod(t(i),t0))<dt/10
        % 迭代
        p(i+1)=p(i)+k*sin(x(i));
    else
        p(i+1)=p(i);
    end      
end

% 绘制粒子运动的动态图像
figure

hold on

% 初始点
h=plot(x(1),0);

% 坐标轴范围
axis([-550 210 -1 1])

% 坐标轴设置
xlabel('x/m');
ylabel('y/m');

% 时刻
title(['t=' num2str(t(1)) 's'])

% 本时刻的位置、速度
text1=text(-20,0.8,['v=' num2str(p(1)/m) 'm/s']);
text2=text(-20,0.6,['x=' num2str(x(1)) 'm']);

for j=2:length(x)
    % 删除上一个点
    delete(h);
    
    % 删除上一个点的位置、速度
    delete(text1);
    delete(text2);
    
    % 绘制下一个点
    h=plot(x(j),0,'ob');
    
    % 时刻
    title(['t=' num2str(t(j)) 's'])
    
    % 本时刻的位置、速度
    text1=text(-20,0.8,['v=' num2str(p(j)/m)  'm/s']);
    text2=text(-20,0.6,['x=' num2str(x(j)) 'm']);
    
    % 动画
    pause(0.1*dt);
end


% 绘制粒子速度、位置随时间的变化
figure

subplot(1,2,1)
plot(t-1e-5,x);

xlabel('t/s');
ylabel('x/m');
title('位置-时间图像');
legend('x-t');


subplot(1,2,2)
plot(t-1e-5,p/m);
xlabel('t/s');
ylabel('v/m\cdot s^{-1}');
title('速度-时间图像');
legend('v-t');


% 绘制相空间轨迹图
figure 
plot(x,p);
xlabel('x/m');
ylabel('p/kg \cdot m/s');
title('动量-时间图像');

