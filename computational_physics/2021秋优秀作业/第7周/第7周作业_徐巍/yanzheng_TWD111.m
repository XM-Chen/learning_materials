%%程序功能：生成某种分布生成的n个随机数，从小到大排列，找到相邻值之差的最大值，
%%看多次重复前面步骤找到的最大值服从的分布是否服从Tracy-Widom分布

clear

n1 = 1000000; %对找1000000个最大值
n2 = 1000; %每次生成的随机数个数
max = zeros(1,n1); %初始化记录最大值的数组

%循环找最大值
for i = 1:n1
    
    %生成各种分布的随机数，需要查看某种分布时，需取消对应的注释符号
    
    %random_number = rand(1,n2); %取值范围为0~1的均匀分布
    random_number = randn(1,n2); %正态分布N(0, 1)
    %random_number = random("Normal",2,5,1, n2); %正态分布N(2, 5)
    %random_number = exprnd(1,1,n2); %指数分布，参数为 1 
    %random_number = gamrnd(1,1,1,n2); %伽马分布，参数为1, 1
    %random_number = raylrnd(1, 1,n2); %瑞利分布，参数为1
    %pd = makedist("rician",1,2); random_number = random(pd, 1, n2); %莱斯分布，参数为1，2
    %pd = makedist("nakagami",1,2); random_number = random(pd, 1, n2); %nakagami分布，参数为1，2
    
    random_number = sort(random_number);%随生成的随机数进行排序
    
    %找出相邻值之差的最大值
    max_dif = 0;
    for j = 1:length(random_number)-1
        if random_number(j+1) - random_number(j) > max_dif
            max_dif = random_number(j+1) - random_number(j);
        end
    end
    max(i) = max_dif;
    
end

%figure()
%计算直方图的总面积，用于归一化
[histFreq, histXout] = hist(max,1000);
binWidth = histXout(2)-histXout(1);
area = binWidth*sum(histFreq);

figure()

%作最大值分布直方图
subplot(2,2,1)
hist(max,1000);
grid on;
title("最大值分布直方图", "FontSize",11)
xlabel('max', "FontSize",11);
ylabel('频数', "FontSize",11);

%作归一化直方图
subplot(2,2,2)
bar(histXout, histFreq/area);
grid on;
title("归一化直方图", "FontSize",11)
xlabel('max', "FontSize",11);
ylabel('频率/组距', "FontSize",11);

%作归一化直方图轮廓图
subplot(2,2,3)
plot(histXout, histFreq/area);
grid on;
title("归一化直方图轮廓", "FontSize",11)
xlabel('max', "FontSize",11);
ylabel('频率/组距', "FontSize",11);

%作平滑轮廓图（近似概率密度分布）
subplot(2,2,4)
%用内核分布平滑轮廓图
PD = fitdist(max', "Kernel");
plot(histXout, pdf(PD, histXout));
grid on
title("平滑轮廓（近似概率密度分布）", "FontSize",11)
xlabel('max', "FontSize",11);
ylabel('频率/组距', "FontSize",11);

ss = 0;
for i = 1: 1000
    ss = binWidth*histFreq(i);
end
