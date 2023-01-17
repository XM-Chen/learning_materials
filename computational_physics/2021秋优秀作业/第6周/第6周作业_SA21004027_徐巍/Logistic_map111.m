%Logistic map数值求解

clear
format long

%x的初值
x(1) = 0.3;

%参数γ
r = 0.025:0.00025:1;

%用迭代式进行迭代
for j = 1:length(r)
    x(1) = 0.3;
    for i = 1:20000
        x(i+1) = 4*r(j)*x(i)*(1-x(i));
    end
    xx(j,:) = x;
end

%作出迭代结果图(稳定值和混沌时的部分值)
figure()
for j = 1:length(r)
    if r(j)<2.9/4
        plot(r(j),xx(j,length(x)),".");%当稳定值只有一个时，只画一个点
        hold on
        continue
    elseif r(j)<3.55/4
        for jj = 19970:20000
            plot(r(j),xx(j,jj),".");
            hold on
            continue
        end
    else
        for jj = 19900:20000
            plot(r(j),xx(j,jj),".");
            hold on
        end
    end
    ylabel('x');
    title("Logistic Map")
    xlabel('r');
end
