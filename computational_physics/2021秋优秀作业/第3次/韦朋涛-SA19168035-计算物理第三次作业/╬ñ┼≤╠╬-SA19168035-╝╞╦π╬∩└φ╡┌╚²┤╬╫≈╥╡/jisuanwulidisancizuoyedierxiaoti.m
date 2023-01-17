% % 计算物理第三次作业第二题
clc;
N = 10001;
L = 10;
h = L/(N-1);
A=0.001;
theta=0.001;
H0 = zeros(N);
V = zeros(N);
for i = 1 : N
    for j = 1 : N
        if i == j
            H0(i,j) = -2;
        else
            if abs(i-j) == 1
                H0(i,j) = 1;
            end
        end
    end
end
x=linspace(-5,5,N);
for i = 1:N
   V(i,i) = ((x(1,i))^2)/2+A*cos(x(1,i)+theta);
end
a=-H0/(2*h*h);
H = -H0/(2*(h^2))+V;
E=eig(H);
% n=linspace(1,N,N);
% plot(n,E)