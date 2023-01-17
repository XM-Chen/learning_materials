clear 
clc

% 实现 Kicked totor model的模拟
% 迭代法
% 改变初始条件k

% 预设
m=1;

% 改变k
n=500;

% 初始条件
x0=1.5;
p0=1;

% 迭代次数
N=5e5;

% 初始化
x=zeros(n,N+1);
p=zeros(n,N+1);

D
for i=1:n
    x(i,1)=mod(x0,2*pi);
    p(i,1)=mod(p0,2*pi);
end

for j=1:n
    k=0.8+0.01*j;

    % 迭代
    for i=1:N
        p(j,i+1)=p(j,i)+k*sin(x(j,i));
        p(j,i+1)=mod(abs(p(j,i+1)),2*pi)*sign(p(j,i+1));
        x(j,i+1)=x(j,i)+p(j,i+1);
        x(j,i+1)=mod(abs(x(j,i+1)),2*pi);
    end
    
end

% 绘制相空间图像
figure
set(gcf,'outerposition',get(0,'screensize'));

hold on

h=plot(x(1,:),p(1,:));

axis([0 2*pi+0.2 -2*pi-0.2 2*pi+0.2])

xlabel('x');
ylabel('p');

title(['k=' num2str(0.11)])


for j=2:n
    delete(h);
    
    h=plot(x(j,:),p(j,:),'.b');
    
    title(['k=' num2str(0.8+0.01*j)])
    
    pause(0.1);
end
