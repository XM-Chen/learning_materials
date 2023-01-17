clear all;
N = 100; %绘图格数
M = 10000; %迭代次数
gamma = 3.83;
xx = linspace(0,1,N+1);
x = zeros(1,M+1);
y = zeros(1,M+1);
y1 = xx;
y2 = gamma*xx.*(1-xx);
plot(xx,y1,'b',xx,y2,'r')
x(1,1) = 0.2; %初值选取
y(1,1) = gamma*x(1,1)*(1-x(1,1));
line([x(1,1),x(1,1)],[0,y(1,1)],'color','g','linewidth',0.1)
for i = 1:M-1
    y(1,i) = gamma*x(1,i)*(1-x(1,i));
    x(1,i+1) = y(1,i);
    line([x(1,i),x(1,i+1)],[y(1,i),x(1,i+1)],'color','g','linewidth',0.1)
    line([x(1,i+1),x(1,i+1)],[x(1,i+1),gamma*x(1,i+1)*(1-x(1,i+1))],'color','g','linewidth',0.1)
end
    