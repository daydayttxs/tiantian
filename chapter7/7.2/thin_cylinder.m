%%%%%%%%%%%�Ŵ��㷨ϡ��Բ����%%%%%%%%%%
%%%%%%%%%%%%%%��ʼ������%%%%%%%%%%%%%%%
clear all;                 %������б���
close all;                 %��ͼ
clc;                       %����
G = 200;                   %����Ŵ�����
NP = 50;                   %��Ⱥ����
Pc = 0.8;                  %������
Pm = 0.05;                 %������
d = 0.5;                   %������Ԫ���,�뱶����
lamda = 1;                 %����
M = 20;                    %����������Բ����
N = 24;                    %ÿ��Բ������Զ��
L = M*N;                   %������Ԫ����
NL = 200;                  %ʵ����Ԫ����              
k = 2;                     %��Ԫ�뾶�벨��֮��      
seta0 = pi/2;              %�źŵ�����λ��           
fai0 = pi;                 %�źŵ���������
NA = 360;                  %�ռ䷽λ�ǲ�����
NE = 360;                  %�ռ丩���ǲ�����
%%%%%%%%%%���ɳ�ʼ��Ⱥ%%%%%%%%%
f = randn(L,NP);
[sortff,Index] = sort(f);
f = ones(L,NP);
for np = 1:NP
    f(Index(end-NL+1:end,np),np) = 0;
end
%%%%%%%%%%%%�Ŵ��㷨ѭ��%%%%%%%%%%%
for g = 1:G
    g
    %%%%%%%%������Ӧ�ȣ�����ֵ�԰�%%%%%%
    for i = 1:NP
        Fit(i) = func_cylinder(k,NE,NA,M,N,seta0,fai0,f(:,i));
    end
    maxFit = max(Fit);           %���ֵ
    minFit = min(Fit);           %��Сֵ
    rr = find(Fit==maxFit);      %�ҳ����ֵ
    fBest = f(:,rr(1,1));        %�������Ÿ���   
    Fit = (Fit-minFit)/(maxFit-minFit);  %��һ����Ӧ��ֵ
    %%%%%%%%%�������̶ĵĸ��Ʋ���%%%%%%%%%
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
    %%%%%%%%%%%���ڸ��ʵĽ������%%%%%%%%
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
    %%%%%%%%%%���ڸ��ʵı������%%%%%%%%
   for m = 1:NP
        for n = 1:L
            r = rand(1,1);
            if r < Pm
                nf(n,m) = ~nf(n,m);
            end
        end
   end
    %%%%%%ʹ��������� ʵ����Ԫ��������%%%%%%
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
    f(:,1) = fBest;                   %�������Ÿ���������Ⱥ��
    trace(g) = maxFit;                %����������Ӧ��
end
figure
plot(trace)
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')
grid on
save fBest.mat fBest