%计算Feigenbaum Constant

clear
format long

%参数γ
r = 0.7:0.0000025:0.95;

j = 1;

%李雅普诺夫指数的绝对值小于该值时记录对应的参数γ
e = 0.0016;

%迭代次数
n = 10000;

%计算各个γ取值下的李雅普诺夫指数，并做判断
for i = 1:length(r)
    s = 0;
    x1 = 0.3;
    
    %用迭代式进行迭代
    for ii = 1:n
        x2 = 4*r(i)*(x1-x1^2);
        s = s + log(4*abs(r(i)-2*r(i)*x1));
        x1 = x2;
    end
    
    %对李雅普诺夫指数做判断
    s = s/n;
    if abs(s) < e 
            rr(j) = r(i);
            ss(j) = s; 
            j = j + 1; 
    end
end

%画出符合要求的各个李雅普诺夫指数
figure()
plot(ss)
xlabel('n');
ylabel('s');


%找到图中前n个最近接0的s值
flag = 1;
nn0 = 5;%只找前五个点，更多的点需要调参数
nn = 1;
jj = 1;
for i = 1:length(ss)
    if nn <= nn0
        if flag ==1
            if ss(i+1)<ss(i)
                sss(jj) = ss(i);
                rrr(jj) = rr(i);
                jj = jj + 1;
                flag = -flag;
                 nn = nn +1;
            end
        else
            if ss(i+1)>ss(i)
                flag = -flag;
            end
        end
    end
end

%在上面找到的5个γ值左右0.0002范围内，分别计算李雅普诺夫指数
j1  =  1;
for i2 = 1:length(rrr)

    clear rr
clear r
clear ss

r = rrr(i2)-0.0002:0.0000005:rrr(i2)+0.0002;

j = 1;
rr(j)= 0;
e = 0.002;
n = 10000;

for i = 1:length(r)
    s = 0;
%     for i1 = 1:50
%     x1 = rand(1);
    x1 = 0.3;
    for ii = 1:n
        x2 = 4*r(i)*(x1-x1^2);
        s = s + log(4*abs(r(i)-2*r(i)*x1));
        x1 = x2;
    end
%     end
    s = s/n;
    if abs(s) < e 
            rr(j) = r(i);
            ss(j) = s; 
            j = j + 1; 
    end
end

%画出符合要求的各个李雅普诺夫指数
figure()
plot(ss)
xlabel('n');
ylabel('s');

%计算前5个分叉点
s0 = 1;
for i = 1:length(ss)
    if abs(ss(i)) < s0
        j = i;
        s0 = abs(ss(i));
    end
end
rrrr(j1) =  rr(j);
j1 = j1+1;
end

%输出前5个分叉点对应的γ值
rrrr

%计算Feigenbaum Constant，并输出
for i = 1:length(rrrr)-2
     (rrrr(i+1)-rrrr(i))/(rrrr(i+2)-rrrr(i+1))
end
