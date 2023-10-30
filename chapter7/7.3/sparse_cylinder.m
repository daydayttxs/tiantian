%%%%%%%%%%%%%遗传算法稀布圆柱阵%%%%%%%%%%%%%%
%%%%%%%%%%%%%%初始化参数%%%%%%%%%%%%%%%
clear all;                %清除所有变量
close all;                %清图
clc;                      %清屏
G = 100;                  %最大遗传代数
NP = 50;                  %种群数量
Pc = 0.8;                 %交叉率
Pm = 0.05;                %变异率
d = 0.5;                  %满阵阵元间距,半倍波长
lamda = 1;                %波长
H=4.5;                    %圆柱阵列高度
R=2;                      %圆柱阵列半径
Ny = 10;                  %俯仰阵元个数
Nx = 12;                  %圆环上阵元个数
L = Nx*Ny;                %稀布阵元个数
k = 2;                    %阵元半径与波长之比 
seta0 = pi/2;             %俯仰指向
fai0 = pi;                %方位指向
NA = 360;                 %空间方位角采样数
NE = 360;                 %空间俯仰角采样数
%%%%%%%%%%生成初始种群%%%%%%%%%
%%%%%%%%%圆环上%%%%%%%%%%
da = acosd((2*R^2-d^2)/(2*R^2));
faimax = 360-Nx*da;
faimin = 0;
Dx = da*ones(Nx,Ny,NP);
Dx(1,:,:) = 0;
Dx = cumsum(Dx,1);
fx = rand(Nx,Ny,NP)*(faimax-faimin)+faimin;   
fx = sort(fx,1);
fx1 = Dx+fx;
%%%%%%%%%俯仰向%%%%%%%%%%
hmax = H-(Ny-1)*d;
hmin = 0;
Dy = d*ones(Nx,Ny,NP);
Dy(:,1,:) = 0;
Dy = cumsum(Dy,2);
fy = rand(Nx,Ny,NP)*(hmax-hmin)+hmin;         
fy = sort(fy,2);
fy1 = Dy+fy;
fy1(1,1,:) = 0;
fy1(1,end,:) = H;
%%%%%%%%%种群%%%%%%%%%%%
f = fy+i.*fx;
ff = fy1+i.*fx1;
for g = 1:G
    g
   %%%%%%计算适应度，即峰值旁瓣%%%%
    for i = 1:NP
        Fit(i) = func_sparse_cylinder(k,NE,NA,seta0,fai0,ff(:,:,i));
    end
    maxFit = max(Fit);           %最大值
    minFit = min(Fit);           %最小值
    rr = find(Fit==maxFit);      %找出最大值
    fBest = ff(:,:,rr(1,1));     %历代最优个体
    Fit = (Fit-minFit)/(maxFit-minFit);  %归一化适应度值
    %%%%%%%基于轮盘赌的复制操作%%%%%%
    sum_Fit = sum(Fit);
    fitvalue = Fit./sum_Fit;
    fitvalue = cumsum(fitvalue);
    ms = sort(rand(NP,1));
    fiti = 1;
    newi = 1;
    while newi<=NP
        if (ms(newi))<fitvalue(fiti)
            nf(:,:,newi) = f(:,:,fiti);
            newi = newi+1;
        else
            fiti = fiti+1;
        end
    end
    %%%%%%%%%%基于概率的交叉操作%%%%%%%%
    for i = 1:2:NP
        p = rand;
        if p<Pc
            q = randi([0,1],1,Nx);
            for j = 1:Nx
                if q(j)==1;
                    temp = nf(j,:,i+1);
                    nf(j,:,i+1) = nf(j,:,i);
                    nf(j,:,i) = temp;
                end
            end
        end
    end
    %%%%%%%基于概率的变异操作%%%%%%%
    for m = 1:NP
        for n = 1:Nx
            r = rand;
            if r<Pm
                nf(n,:,m) = rand(1,Ny)*(hmax-hmin)+hmin...
                    +sqrt(-1)*(rand(1,Ny)*(faimax-faimin)+faimin);
            end
        end
    end
    fx = imag(nf);    
    fx = sort(fx,1);
    fx1 = Dx+fx; 
    fy = real(nf);
    fy = sort(fy,2);
    fy1 = Dy+fy; 
    fy1(1,1,:) = 0;
    fy1(1,end,:) = H; 
    f = fy+sqrt(-1).*fx;
    ff = fy1+sqrt(-1).*fx1;
    ff(:,:,1) = fBest;                %保留最优个体在新种群中
    trace(g) = maxFit;                %历代最优适应度
end
figure
plot(trace)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')
grid on
save fBest.mat fBest