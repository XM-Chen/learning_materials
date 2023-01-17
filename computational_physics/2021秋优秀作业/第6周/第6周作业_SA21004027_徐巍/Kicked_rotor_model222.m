%数值求解Kicked rotor model
clear

K = 1.5;%参数k
n=10%取10组初始值，并将结果作在一张图中


for jj = 1:n
    
N = 20000;%迭代次数

%随机取初值
x1 =  rand(1)*2*pi;
p1 =  rand(1)*2*pi;

%数值计算
for i = 1:N
    x(i) = x1;
    p(i) = p1;
   p1 =p1+K*sin(x1);
    x1 = x1+p1;
end

pp(jj,:)=p;
xx(jj,:)=x;
end

%在相空间可视化10组初始值的求解结果
figure()
for i = 1:n
plot(mod(xx(i,:),2*pi),mod(pp(i,:),2*pi),".")
hold on
xlabel('x');
ylabel('p');
title("k = "+num2str(K)+"（相空间）")
axis([0,2*pi,0,2*pi])
end
