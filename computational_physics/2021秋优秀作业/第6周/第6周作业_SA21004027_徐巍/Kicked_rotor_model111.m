%数值求解Kicked rotor model

clear

K = 0.1;%参数k
n = 1%取1组初始值


%迭代次数
N = 20000;

%初值
x1 =  1;
p1 =  1;

%数值计算
for i = 1:N
    x(i) = x1;
    p(i) = p1;
   p1 =p1+K*sin(x1);
   x1 = x1+p1;
end

%在普通空间可视化求解结果
figure()
plot(x,p,".")
xlabel('x');
ylabel('p');
title("k = "+num2str(K)+"（普通空间）")
hold on

%在相空间可视化求解结果
figure()
for i = 1:n
plot(mod(x,2*pi),mod(p,2*pi),".")
hold on
xlabel('x');
ylabel('p');
title("k = "+num2str(K)+"（相空间）")
axis([0,2*pi,0,2*pi])
end



