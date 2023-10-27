%%�Ŵ��㷨��ϡ������
%%%%ȫ������Ԫ128����ȡ64��active��Ԫ��ϡ������
clear all;
close all;
clc;
c=1540;
f=3e6;
NP =112;  %��Ⱥ����
Pc = 0.8;   %������
Pm = 0.03;   %������
lamda = c/f;     %����
 d = 0.9*lamda;      %������Ԫ���
% d = 0.108e-3;      %������Ԫ���
seta0 = 0*pi/180;   %����ָ��
G = 100;         %����Ŵ�����   ����������
L = 64;          %������Ԫ����
NL = 32;         %ʵ����Ԫ����
NN = 1800;        %ɨ��Ƕȼ��  ������Ԫ��ȷ��
%%%%%%%%%%%%%%%���ɳ�ʼ��Ⱥ%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%1Ϊ����Ԫ��0Ϊ����Ԫ%%%%%%%%%%%%%%%%%
f = randn(L,NP);
[sortff,Index] = sort(f);
f = zeros(L,NP);
for i=1:NP
    f(Index(end-NL+1:end,i),i)=1;
end
%%%%%%%%%%%%%%% �Ŵ��㷨ѭ�� %%%%%%%%%%%%%%%
for k=1:G
    k
    %%%%%%%%%%%%%%% ������Ӧ�ȣ�����ֵ�԰�� %%%%%%%%%%%%%%%
    for i=1:NP
        Fit(i) = -func_line(d,lamda,NN,seta0,f(:,i));
    end
    maxFit = max(Fit);    %���ֵ
    minFit = min(Fit);      %��Сֵ
    rr = find(Fit==maxFit);    %�ҳ����ֵ
    fBest = f(:,rr(1,1));      %�������Ÿ���
    Fit = (Fit-minFit)/(maxFit-minFit);    %��һ����Ӧ��ֵ
 %%%%%%%%%%%%%%% �������̶ĵ�ѡ����� %%%%%%%%%%%%%%%
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
 %%%%%%%%%%%% ���ڸ��ʵĽ������%%%%%%%%%
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
 %%%%%%%%%%%% ���ڸ��ʵı������%%%%%%%%%
 for m = 1: NP
     for n = 1:L
            r = rand(1,1);
            if r < Pm
                nf(n,m)= ~nf(n,m);
            end
     end
 end
 %%%%%%%%%%%% ʹ���������ʵ����Ԫ��������%%%%%%%%%
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
  f(:,1) = fBest;     %�������Ÿ���������Ⱥ��
  trace(k) = maxFit;     %����������Ӧ��ֵ
end

plot(trace,'r')
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������/��ֵ�԰��')
grid on 
save fBest3.mat fBest
