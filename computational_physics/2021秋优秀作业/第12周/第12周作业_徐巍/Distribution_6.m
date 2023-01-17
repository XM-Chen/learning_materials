%%第（3）题情况（b)
clear

N = 1000;
n11 = 10000;%重复计算次数
for i1 = 1:n11
i1
M = randn(N);

for ii = 1:length(M)
    for jj = 1:length(M)
        if jj > ii
            M(ii, jj) = normrnd(0, 1, 1);
            M(jj, ii) = conj(M(ii, jj));
        end
    end
end


[x, y] = eig(M);

eigenvalues = diag(y);%已经按从小到大排序

for i = 1:length(eigenvalues) - 1;
    dd(i) = eigenvalues(i+1) - eigenvalues(i);
end

maxD(i1) = max(dd);
end


nbins = 100;%直方图的区间数


[histFreq, histXout]=hist(maxD, nbins);

binWidth = histXout(2)-histXout(1);
area = binWidth*sum(histFreq);

figure()
subplot(3,1,1)
bar(histXout,histFreq/area);
grid on;
title("归一化直方图")
xlabel("Δ")
ylabel("P")


subplot(3,1,2)
plot(histXout,histFreq/area,'r*');
grid on;
title("直方图顶点散点图")
xlabel("Δ")
ylabel("P")
 
subplot(3,1,3)
bar(histXout,histFreq/area);
hold on;
plot(histXout,histFreq/area,'r-',"LineWidth",2);
grid on;
title("直方图和散点图")
xlabel("Δ")
ylabel("P")

x = histXout;
y = histFreq/area;

[xData, yData] = prepareCurveData(x, y);

%按指定函数形式拟合
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [2.17292572042206, 0.846164873822668, 0.0966154968178991, 1.38253167469357, 0.969869056378798, 0.1242249338958];


[fitresult, gof] = fit( xData, yData, ft, opts );
fitresult

a1 = fitresult.a1;
b1 = fitresult.b1;
c1 = fitresult.c1;
a2 = fitresult.a2;
b2 = fitresult.b2;
c2 = fitresult.c2;


xmin = min(histXout);
xmax = max(histXout);

x1 = xmin:0.001:xmax;
y1 = a1*exp(-((x1-b1)/c1).^2) + a2*exp(-((x1-b2)/c2).^2);
figure()
plot(x1,y1,"b-","LineWidth",1.5)
hold on
plot(x,y,"r.","MarkerSize",11)
legend("拟合曲线", "数值散点")
xlabel("Δ")
ylabel("P(Δ)")
