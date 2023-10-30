%%%%%%%%%%%遗传算法稀疏圆柱阵%%%%%%%%%%
%%%%%%%%%%%%%%初始化参数%%%%%%%%%%%%%%%
clear all;                 %清除所有变量
close all;                 %清图
clc;                       %清屏
G = 200;                   %最大遗传代数
NP = 50;                   %种群数量
Pc = 0.8;                  %交叉率
Pm = 0.05;                 %变异率
d = 0.5;                   %满阵阵元间距,半倍波长
lamda = 1;                 %波长
M = 20;                    %柱面阵列中圆弧数
N = 24;                    %每个圆弧中阵远数
L = M*N;                   %满阵阵元个数
NL = 200;                  %实际阵元个数              
k = 2;                     %阵元半径与波长之比      
seta0 = pi/2;              %信号到来方位角           
fai0 = pi;                 %信号到来俯仰角
NA = 360;                  %空间方位角采样数
NE = 360;                  %空间俯仰角采样数
%%%%%%%%%%生成初始种群%%%%%%%%%
f = randn(L,NP);
[sortff,Index] = sort(f);
f = ones(L,NP);
for np = 1:NP
    f(Index(end-NL+1:end,np),np) = 0;
end
%%%%%%%%%%%%遗传算法循环%%%%%%%%%%%
for g = 1:G
    g
    %%%%%%%%计算适应度，即峰值旁瓣%%%%%%
    for i = 1:NP
        Fit(i) = func_cylinder(k,NE,NA,M,N,seta0,fai0,f(:,i));
    end
    maxFit = max(Fit);           %最大值
    minFit = min(Fit);           %最小值
    rr = find(Fit==maxFit);      %找出最大值
    fBest = f(:,rr(1,1));        %历代最优个体   
    Fit = (Fit-minFit)/(maxFit-minFit);  %归一化适应度值
    %%%%%%%%%基于轮盘赌的复制操作%%%%%%%%%
    sum_Fit = sum(Fit);
    fitvalue = Fit./sum_Fit;
    fitvalue = cumsum(fitvalue);
    ms = sort(rand(NP,1));
    fiti = 1;
    newi = 1;
    while newi<=NP
        if (ms(newi))<fitvalue(fiti)
            nf(:,newi) = f(:,fiti);
            newi = newi+1;
        else
            fiti = fiti+1;
        end
    end   
    %%%%%%%%%%%基于概率的交叉操作%%%%%%%%
    for i = 1:2:NP
        p = rand;
        if p<Pc
            q = randi([0,1],1,L);
            for j = 1:L
                if q(j)==1;
                    temp = nf(j,i+1);
                    nf(j,i+1) = nf(j,i);
                    nf(j,i) = temp;
                end
            end
        end
    end
    %%%%%%%%%%基于概率的变异操作%%%%%%%%
   for m = 1:NP
        for n = 1:L
            r = rand(1,1);
            if r < Pm
                nf(n,m) = ~nf(n,m);
            end
        end
   end
    %%%%%%使交叉变异后的 实际阵元个数不变%%%%%%
    for i = 1:NP
        n_ones = sum(nf(:,i));
        while n_ones>(NL)   
            nn1 = find(nf(:,i)==1);
            MUT1 = randi([1,n_ones],1,n_ones-NL);
            for m = 1:(n_ones-NL)
                nf(nn1(MUT1(m)),i) = 0;
            end
            n_ones = sum(nf(:,i));
        end
        while n_ones<(NL)
            nn2 = find(nf(:,i)==0);
            MUT2 = randi([1,L-n_ones],1,NL-n_ones);
            for m = 1:(NL-n_ones)
                nf(nn2(MUT2(m)),i) = 1;
            end
            n_ones = sum(nf(:,i));
        end
    end
    f = nf;
    f(:,1) = fBest;                   %保留最优个体在新种群中
    trace(g) = maxFit;                %历代最优适应度
end
figure
plot(trace)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')
grid on
save fBest.mat fBest