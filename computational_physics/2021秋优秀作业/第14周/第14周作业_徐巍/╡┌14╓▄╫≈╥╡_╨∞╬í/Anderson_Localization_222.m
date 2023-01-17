%%Anderson Localization——题目第二问，V*cos(arfa*n1*pi)
clear
format long

t = 1;

n = 3000;

V0 = [0,0.4,0.7,0.8,0.9,1,1.1,1.2,1.3,2,3,4,5]*2;
%V0 = 0:0.2:30

eg = [];
for ii = 1:length(V0)
V = V0(ii);
arfa = (sqrt(5)-1)/2;
n1 = 1:n;
h1 = V*cos(arfa*n1*pi);


h2 = ones(1,n-1)*(-t);

H1 = diag(h1);
H2 = diag(h2,1);
H3 = diag(h2,-1);
H = H1 + H2 + H3;




[x, y] = eig(H);
eigenvalues = diag(y);
eg = [eg,eigenvalues];


figure()
 subplot(2,1,1)
 plot(x(:,1));
 hold on
 plot(x(:,1500));
 plot(x(:,2998));
 xlabel("Lattice sites")
 ylabel("phi")
 title("波函数"+"(V = "+num2str(V0(ii))+")")
 legend("基态", "第1500个本征态","第2998个本征态")
 
 subplot(2,1,2)
 plot(x(:,1).^2);
 hold on
 plot(x(:,1500).^2);
 plot(x(:,2998).^2);
 xlabel("Lattice sites")
 ylabel("|phi|^2")
 title("态密度图"+"(V = "+num2str(V0(ii))+")")
 legend("基态", "第1500个本征态","第2998个本征态")

 figure()
 plot(eigenvalues,".")
 xlabel("第n个态")
 ylabel("En")
 title("能谱"+"(V = "+num2str(V0(ii))+")")

end


clear x
figure()

for i = 1:length(V0)
    x = ones(length(eigenvalues),1)*V0(i);
    plot(x,eg(:,i),".")
    hold on
end
xlabel("V")
ylabel("En")
title("能谱")

