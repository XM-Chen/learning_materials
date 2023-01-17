clear 
clc

% 用于求解logistic map问题
% 使用函数迭代


% ni是分割lambda的格点数
% nj为函数迭代次数
ni=200;
nj=700;

% 对x进行初始化
% 列代表不同的lambda
% 行代表迭代次数
x=zeros(nj,ni);
for i=1:ni
    % lambda的大小
    lambda=0.72+0.001*(i-1);
    
    % 函数定义
    f=@(x) 4*lambda.*x.*(1-x);
    
    % 进行函数迭代
    % 从x=0.5开始迭代
    for j=1:nj
        if j==1
            x(j,i)=f(0.5);
        else
            x(j,i)=f(x(j-1,i));
        end
    end
end

% 将lambda的数值转化为矩阵形式
lambda=0.72:0.001:0.919;

% 进行绘图
figure
hold on

% 图例设置
xlabel('\lambda');
ylabel('x');
title('logistic map');

% 绘制散点图
for i=600:nj
    
    % 只绘制后100个点
    % 迭代次数较多
    % 稳定点已稳定
    % 分叉点已出现
    
    scatter(lambda,x(i,:))
end


