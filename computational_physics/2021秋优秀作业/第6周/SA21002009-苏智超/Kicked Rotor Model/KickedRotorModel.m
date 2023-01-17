clear all;
k = 0.7;
N = 1000000;
x = zeros(1,N);
p = zeros(1,N);
x(1,1) = 0.5;
p(1,1) = 0.1;
for i = 2:N
    p(1,i) = p(1,i-1)+k*sin(x(1,i-1));
    x(1,i) = x(1,i-1)+p(1,i);
end
plot(x,p,'.')
xlabel('x');
ylabel('p');