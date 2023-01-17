%%Anderson Localization——题目第三问，V*cos(arfa*abs(n1).^v)
clear
format long

t = 1;

n = 3000;

% V0 = [0,0.7,1,1.3,2]*2;
% arfa = 1;
% v = 0.3;

% V = 1;
% arfa0 = 0.1:0.1:0.5;
% v = 0.3;

V = 1;
arfa = 0.3;
v0 = 0.1:0.1:0.5;

eg = [];

for ii = 1:length(v0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = v0(ii)

n1 = 1:n;
h1 = V*cos(arfa*abs(n1).^v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
 title("波函数"+"(v = "+num2str(v0(ii))+")")
 legend("基态", "第1500个本征态","第2998个本征态")
 
 subplot(2,1,2)
 plot(x(:,1).^2);
 hold on
 plot(x(:,1500).^2);
 plot(x(:,2998).^2);
 xlabel("Lattice sites")
 ylabel("|phi|^2")
 title("态密度图"+"(v = "+num2str(v0(ii))+")")
 legend("基态", "第1500个本征态","第2998个本征态")

 figure()
 plot(eigenvalues,".")
 xlabel("第n个态")
 ylabel("En")
 title("能谱"+"(v = "+num2str(v0(ii))+")")

end


clear x
figure()

for i = 1:length(v0)
    x = ones(length(eigenvalues),1)*v0(i);
    plot(x,eg(:,i),".")
    hold on
end
xlabel("V")
ylabel("En")
title("能谱")
