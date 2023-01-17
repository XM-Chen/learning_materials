% % 本程序为中科大2021年秋季学期计算物理第三次作业第三小题(b)
% % 用数值方法求解心形无限深势阱的本征值、本征态
% % 
clear;clc;
v=20000;
h=0.1;

x=linspace(-2.5,2.5,51);
y=linspace(-3.5,1.5,51);
x_coordinate=diag(kron(eye(length(y)),diag(x)));% 每一行的横坐标，每行51个
y_coordinate=diag(kron(diag(y),eye(length(x))));% 每一行的纵坐标，每行51个

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
    if (y0<=sqrt(2*sqrt(x0^2)-x0^2))&&(y0>=asin(abs(x0)-1)-pi/2)
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

ev=EV(:,4);% 获取本征矢

xx=-2:0.01:2;
yy=sqrt(2*sqrt(xx.^2)-xx.^2);% 心形的上部分
zz=asin(abs(xx)-1)-pi./2;% 心形的下部分
plot(xx,yy,'red')
grid on
hold on
plot(xx,zz,'red')
axis([-3,3,-4,2])
hold on

theta=0:2*pi;
for i = 1 : N
%     scatter(x_coordinate(i,1),y_coordinate(i,1),'.','black')% 画散点图
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
