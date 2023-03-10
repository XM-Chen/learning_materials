% % 本程序为中科大2021年秋季学期计算物理第三次作业第三小题(a)
% % 用数值方法求解椭圆形无限深势阱的本征值、本征态
% % 
clear;clc;
h=0.10000;
v=20000;
a=5;% 长轴
b=3;% 短轴

x=linspace(-5.1,5.1,103);
y=linspace(-3.1,3.1,63);
x_coordinate=diag(kron(eye(length(y)),diag(x)));% 每一行的横坐标，每行103个
y_coordinate=diag(kron(diag(y),eye(length(x))));% 每一行的纵坐标，每行103个

N=length(x)*length(y);

V=zeros(N);
H0=zeros(N);

for i=1:N
    for j=1:N
        if i==j
           H0(i,j) =-4;
        else
            x0=x_coordinate(i,1);
            y0=y_coordinate(i,1);
            x1=x_coordinate(j,1);
            y1=y_coordinate(j,1);
            delta_x=x1-x0;
            delta_y=y1-y0;
            if (sqrt(delta_x^2+delta_y^2)>0.09) && (sqrt(delta_x^2+delta_y^2)<0.11)
                H0(i,j)=1;
            else
                H0(i,j)=0;
            end
        end
    end
end

for i=1:N
    x0=x_coordinate(i,1);
    y0=y_coordinate(i,1);
    if ((x0/a)^2+(y0/b)^2)<=1
        V(i,i)=0;
    else
        V(i,i)=v;
    end
end

H=-H0/(2*h^2)+V;
[EV,E]=eig(H);
E0=diag(E);
index=linspace(1,N,N);
% scatter(index,diag(E))% 画能谱
% xlabel('state index')
% ylabel('Energy')

ev=EV(:,1);% 获取本征矢
% 画椭圆
beta=linspace(0,2*pi,1000);
m=3*sin(beta);
n=5*cos(beta);
plot(n,m,'red')
hold on;

% 画波函数的分布
theta=0:2*pi;
for i=1:N
%     scatter(x_coordinate(i,1),y_coordinate(i,1),'.','black')
    d=abs(ev(i,1));
    if ev(i,1)>0
        x=x_coordinate(i,1)+d*cos(theta);
        y=y_coordinate(i,1)+d*sin(theta);
        fill(x,y,'red')
        plot(x,y,'red')
    else
        x=x_coordinate(i,1)+d*cos(theta);
        y=y_coordinate(i,1)+d*sin(theta);
        fill(x,y,'black')
        plot(x,y,'black')
    end
    if i<=N
        hold on;
    end
end