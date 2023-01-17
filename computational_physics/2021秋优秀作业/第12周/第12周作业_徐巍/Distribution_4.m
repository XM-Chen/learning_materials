%%三维复厄米随机矩阵特征值分布
clear
nnn = 100000000;
for ii = 1:nnn
ii
a11 = randn(1);
a22 = randn(1);
a33 = randn(1);

a12 = normrnd(0,0.5,1) + i*normrnd(0,0.5,1);
a13 = normrnd(0,0.5,1) + i*normrnd(0,0.5,1);
a23 = normrnd(0,0.5,1) + i*normrnd(0,0.5,1);

a21 = conj(a12);
a31 = conj(a13);
a32 = conj(a23);


H = [a11, a12, a13; a21, a22, a23; a31, a32, a33];

[x,y] = eig(H);


r = zeros(1,3);
r(1) = y(1,1);
r(2) = y(2,2);
r(3) = y(3,3);

r = r(randperm(length(r)));


r1(ii) = r(1);
r2(ii) = r(2);
r3(ii) = r(3);

end

figure()
plot(r1,".")
hold on
plot(r2,".")
plot(r3,".")
legend("特征值1", "特征值2","特征值3")
xlabel("循环次数");
ylabel("特征值")


n1 = 6;
bc = 0.1;
rr = -n1:bc:n1;
nr = length(rr);
Ds = zeros(nr,nr,nr);

n0 = 0;
%统计特征值的分布情况
for i = 1:nnn
    
     if r1(i)>6||r1(i)<-6||r2(i)>6||r2(i)<-6||r3(i)>6||r3(i)<-6
        n0 = 1 + n0;
        continue;
    end
    
    ii =  round((r1(i)+n1)*1/bc)+1;
    jj =  round((r2(i)+n1)*1/bc)+1;
    kk =  round((r3(i)+n1)*1/bc)+1;
    
    Ds(ii, jj, kk) = Ds(ii, jj, kk) + 1;
    
end

Ds1 = Ds;

Ds = Ds/(sum(sum(sum(Ds)))*bc^3);

X = rr;
Y = X;
Z = X;


qp = 76;


%%%%固定x时的切片
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xi = qp;
rr(xi)
n10 = length(Ds);
yz =  zeros(n10);
for ii = 1:n10
    for jj = 1:n10
        yz(ii,jj) = Ds(xi, ii, jj);
    end
end


[y, z] = meshgrid(rr,rr);

x1 = y;
y1 = z;
[xData, yData, zData] = prepareSurfaceData( x1, y1, yz );

%按指定函数形式拟合
%ft中的1.5即为rr(xi)的值，每次修改rr(xi)后，需手动改变该值
ft = fittype( 'C*(abs(x-y)*abs(x-(1.5))*abs(y-(1.5)))^w*exp(-A*(x^2+y^2+(1.5)^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.588104609027275, 0.0203237944133983, 0.9080361421463];

[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
fitresult


Const = fitresult.C;
A = fitresult.A;
w = fitresult.w;

P = Const*(abs(rr(xi)-x1).*abs(rr(xi)-y1).*abs(y1-x1)).^w.*exp(-A*(rr(xi)^2 + y1.^2 + x1.^2));

chazhi = yz - P;



figure()
plot3(xData, yData, zData ,"r.","MarkerSize",3)
xlabel("x")
ylabel("y")
zlabel("P")
title("数值统计结果(固定x)")

figure()
mesh(Y, Z, P)
xlabel("y")
ylabel("z")
zlabel("P")
title("分布函数拟合结果(固定x)")


figure()
fig1 = plot3(xData, yData, zData ,"r.","MarkerSize",11);
hold on
fig2 = mesh(Y, Z, P);
xlabel("y")
ylabel("z")
zlabel("P")
title("数值统计结果与分布函数拟合结果对比(固定x)")
legend([fig1, fig2],"数值统计结果","分布函数拟合结果")


figure()
mesh(Y,Z,chazhi)
xlabel("y")
ylabel("z")
zlabel("dP")
title("数值统计结果与分布函数拟合结果之差(固定x)")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%固定y时的切片
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yi = qp;
rr(yi)
n10 = length(Ds);
xz =  zeros(n10);
for ii = 1:n10
    for jj = 1:n10
        xz(ii,jj) = Ds(ii, yi, jj);
    end
end

[x, z] = meshgrid(rr,rr);

x1 = x;
y1 = z;
[xData, yData, zData] = prepareSurfaceData( x1, y1, xz );

% Set up fittype and options.
ft = fittype( 'C*(abs(x-y)*abs(x-(1.5))*abs(y-(1.5)))^w*exp(-A*(x^2+y^2+(1.5)^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.588104609027275, 0.0203237944133983, 0.9080361421463];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
fitresult

Const = fitresult.C;
A = fitresult.A;
w = fitresult.w;

P = Const*(abs(rr(yi)-x1).*abs(rr(yi)-y1).*abs(y1-x1)).^w.*exp(-A*(rr(yi)^2 + y1.^2 + x1.^2));

chazhi = xz - P;



figure()
plot3(xData, yData, zData ,"r.","MarkerSize",3)
xlabel("x")
ylabel("z")
zlabel("P")
title("数值统计结果(固定y)")

figure()
mesh(X, Z, P)
xlabel("x")
ylabel("z")
zlabel("P")
title("分布函数拟合结果(固定y)")


figure()
fig1 = plot3(xData, yData, zData ,"r.","MarkerSize",11);
hold on
fig2 = mesh(X, Z, P);
xlabel("x")
ylabel("z")
zlabel("P")
title("数值统计结果与分布函数拟合结果对比(固定y)")
legend([fig1, fig2],"数值统计结果","分布函数拟合结果")


figure()
mesh(X,Z,chazhi)
xlabel("x")
ylabel("z")
zlabel("dP")
title("数值统计结果与分布函数拟合结果之差(固定y)")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%固定z时的切片
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zi = qp;
rr(zi)
n10 = length(Ds);
xy =  zeros(n10);
for ii = 1:n10
    for jj = 1:n10
        xy(ii,jj) = Ds(ii, jj, zi);
    end
end

[x, y] = meshgrid(rr,rr);

x1 = x;
y1 = y;
[xData, yData, zData] = prepareSurfaceData( x1, y1, xy );

ft = fittype( 'C*(abs(x-y)*abs(x-(1.5))*abs(y-(1.5)))^w*exp(-A*(x^2+y^2+(1.5)^2))', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.588104609027275, 0.0203237944133983, 0.9080361421463];

[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
fitresult

Const = fitresult.C;
A = fitresult.A;
w = fitresult.w;

P = Const*(abs(rr(zi)-x1).*abs(rr(zi)-y1).*abs(y1-x1)).^w.*exp(-A*(rr(zi)^2 + y1.^2 + x1.^2));

chazhi = xy - P;



figure()
plot3(xData, yData, zData ,"r.","MarkerSize",3)
xlabel("x")
ylabel("y")
zlabel("P")
title("数值统计结果(固定z)")

figure()
mesh(X, Y, P)
xlabel("x")
ylabel("y")
zlabel("P")
title("分布函数拟合结果(固定z)")


figure()
fig1 = plot3(xData, yData, zData ,"r.","MarkerSize",11);
hold on
fig2 = mesh(X, Y, P);
xlabel("x")
ylabel("y")
zlabel("P")
title("数值统计结果与分布函数拟合结果对比(固定z)")
legend([fig1, fig2],"数值统计结果","分布函数拟合结果")


figure()
mesh(X,Y,chazhi)
xlabel("x")
ylabel("y")
zlabel("dP")
title("数值统计结果与分布函数拟合结果之差(固定z)")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    
    