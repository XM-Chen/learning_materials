clear all;
M = 1000; %迭代次数
gamma_min = 3.4;
gamma_max = 4;
N = 1000; %绘图格数
num = 500; %初值选取次数
gamma = linspace(gamma_min,gamma_max,N+1);
x = zeros(1,N+1);
a = zeros(1,M);
for k = 1:num-1
    a(1,1) = 1/num*k;
    for i = 1:N+1
        for j = 1:M-1
            a(1,j+1) = gamma(1,i)*a(1,j)*(1-a(1,j));
        end
        x(1,i) = a(1,M);
    end
    plot(gamma,x,'.r','MarkerSize',0.01)
    hold on
end
xlabel('gamma');
ylabel('x');