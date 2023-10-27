%%遗传算法，稀疏阵列
%%%%全阵列阵元128个，取64个active阵元的稀疏阵列
clear all;
close all;
clc;
c=1540;
f=3e6;
NP =112;  %种群数量
Pc = 0.8;   %交叉率
Pm = 0.03;   %变异率
lamda = c/f;     %波长
 d = 0.9*lamda;      %满阵阵元间距
% d = 0.108e-3;      %满阵阵元间距
seta0 = 0*pi/180;   %波束指向
G = 100;         %最大遗传代数   最大迭代次数
L = 64;          %满阵阵元个数
NL = 32;         %实际阵元个数
NN = 1800;        %扫描角度间隔  根据阵元数确定
%%%%%%%%%%%%%%%生成初始种群%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%1为有阵元，0为无阵元%%%%%%%%%%%%%%%%%
f = randn(L,NP);
[sortff,Index] = sort(f);
f = zeros(L,NP);
for i=1:NP
    f(Index(end-NL+1:end,i),i)=1;
end
%%%%%%%%%%%%%%% 遗传算法循环 %%%%%%%%%%%%%%%
for k=1:G
    k
    %%%%%%%%%%%%%%% 计算适应度，即峰值旁瓣比 %%%%%%%%%%%%%%%
    for i=1:NP
        Fit(i) = -func_line(d,lamda,NN,seta0,f(:,i));
    end
    maxFit = max(Fit);    %最大值
    minFit = min(Fit);      %最小值
    rr = find(Fit==maxFit);    %找出最大值
    fBest = f(:,rr(1,1));      %历代最优个体
    Fit = (Fit-minFit)/(maxFit-minFit);    %归一化适应度值
 %%%%%%%%%%%%%%% 基于轮盘赌的选择操作 %%%%%%%%%%%%%%%
    sum_Fit = sum(Fit);
    fitvalue = Fit./sum_Fit;
    fitvalue = cumsum(fitvalue);
    ms = sort(rand(NP,1));
    fiti = 1;
    newi = 1;
    while newi <= NP
        if(ms(newi)) < fitvalue(fiti)
            nf(:,newi) = f(:,fiti);
            newi = newi +1;
        else
            fiti=fiti+1;
        end
     end
 %%%%%%%%%%%% 基于概率的交叉操作%%%%%%%%%
 for i=1:2:NP
     p = rand;
     if p < Pc
         q = randi([0,1],1,L);
         for j= 1:L
             if q(j)==1
                 temp = nf(j, i+1);
                 nf(j, i+1) =  nf(j, i);
                 nf(j, i)=temp;
             end
         end
     end
 end
 %%%%%%%%%%%% 基于概率的变异操作%%%%%%%%%
 for m = 1: NP
     for n = 1:L
            r = rand(1,1);
            if r < Pm
                nf(n,m)= ~nf(n,m);
            end
     end
 end
 %%%%%%%%%%%% 使交叉变异后的实际阵元个数不变%%%%%%%%%
  for i = 1:NP
      n_ones = sum(nf(:,i));
      while n_ones>(NL)
          nn1=find(nf(:,i)==1);
          MUT1 = randi([1,n_ones],1,n_ones-NL);
          for m =1:(n_ones-NL)
              nf(nn1(MUT1(m)),i)=0;
          end
          n_ones = sum(nf(:,i));
      end
      while n_ones<(NL)
          nn2 = find(nf(:,i)==0);
          MUT2= randi([1,L-n_ones],1,NL-n_ones);
          for m = 1:(NL-n_ones)
             nf(nn2(MUT2(m)),i)=1; 
          end
          n_ones = sum(nf(:,i));
      end
  end
  f=nf;
  f(:,1) = fBest;     %保留最优个体在新种群中
  trace(k) = maxFit;     %历代最优适应度值
end

plot(trace,'r')
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线/峰值旁瓣比')
grid on 
save fBest3.mat fBest
