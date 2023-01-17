clear 
clc

% 实现 Kicked totor model的模拟
% 迭代法
% 改变初始条件x0

% 预设
k=1;
m=1;

% 改变x0
n=100;

% 初始条件
x0=2:0.1:2+0.1*(n-1);
p0=1;

% 迭代次数
N=1e5;

% 初始化
x=zeros(n,N+1);
p=zeros(n,N+1);

for i=1:n
    x(i,1)=mod(abs(x0(i)),2*pi)*sign(x0(i));
    p(i,1)=mod(abs(p0),2*pi)*sign(p0);
end

for j=1:n

    % 迭代
    for i=1:N
        p(j,i+1)=p(j,i)+k*sin(x(j,i));
        p(j,i+1)=mod(abs(p(j,i+1)),2*pi)*sign(p(j,i+1));
        x(j,i+1)=x(j,i)+p(j,i+1);
        x(j,i+1)=mod(abs(x(j,i+1)),2*pi);
    end
    
end

% for j=1:n
%     for i=1:N+1
%         x(j,i)=mod(x(j,i),2*pi);
%     end
% end


% 绘制相空间图像
figure
set(gcf,'outerposition',get(0,'screensize'));

hold on

h=scatter(x(1,:),p(1,:));

axis([-0.2 2*pi+0.2 -4 4])

xlabel('x');
ylabel('p');

title(['x0=' num2str(x0(1))])


for j=2:n
    delete(h);
    
    h=scatter(x(j,:),p(j,:),'.b');
    
    title(['x0=' num2str(x0(j))])
    
    pause(0.1);
end
