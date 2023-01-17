% 本程序随机生成Google矩阵并求对应的本征值。
clear;clc;
N=3000; %维数
alpha=0.1; %阻尼系数
mean=10; %平均每个网页的连接数
M_matrix=zeros(N,N); %Markov矩阵
for i=1:N
    for j=1:N
        if(rand<mean/N)
            M_matrix(i,j)=1;
        end
    end
end
for i=1:N
    if(sum(M_matrix(i,:))~=0)
        M_matrix(i,:)=M_matrix(i,:)/sum(M_matrix(i,:));
    else
        M_matrix(i,:)=1/N*ones(1,N);
    end
end
G_matrix=(alpha*M_matrix+(1-alpha)/N*ones(N,N))';
[vector,value]=eig(G_matrix);
[spectrum,index]=sort(diag(value),'descend');
vector_sort=vector(:,index);
probability=abs(vector(:,1));
figure(1)
plot(1:N,abs(spectrum),'r');
title('The module of the eigenvalues of Google matrix(alpha=0.1)');
xlabel('label $i$','interpreter','latex');
ylabel('$|\lambda|_i$','interpreter','latex');
figure(2)
scatter(real(spectrum),imag(spectrum),2,'b','filled');
title('The spectrum of the eigenvalues of Google matrix(alpha=0.1)');
xlabel('Re($\lambda$)','interpreter','latex');
ylabel('Im($\lambda$)','interpreter','latex');
page_rank=abs(sort(vector_sort(:,1),'descend'));
% figure(3)
% plot((1:2000),log10(page_rank),'g');