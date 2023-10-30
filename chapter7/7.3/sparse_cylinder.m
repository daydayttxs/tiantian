%%%%%%%%%%%%%�Ŵ��㷨ϡ��Բ����%%%%%%%%%%%%%%
%%%%%%%%%%%%%%��ʼ������%%%%%%%%%%%%%%%
clear all;                %������б���
close all;                %��ͼ
clc;                      %����
G = 100;                  %����Ŵ�����
NP = 50;                  %��Ⱥ����
Pc = 0.8;                 %������
Pm = 0.05;                %������
d = 0.5;                  %������Ԫ���,�뱶����
lamda = 1;                %����
H=4.5;                    %Բ�����и߶�
R=2;                      %Բ�����а뾶
Ny = 10;                  %������Ԫ����
Nx = 12;                  %Բ������Ԫ����
L = Nx*Ny;                %ϡ����Ԫ����
k = 2;                    %��Ԫ�뾶�벨��֮�� 
seta0 = pi/2;             %����ָ��
fai0 = pi;                %��λָ��
NA = 360;                 %�ռ䷽λ�ǲ�����
NE = 360;                 %�ռ丩���ǲ�����
%%%%%%%%%%���ɳ�ʼ��Ⱥ%%%%%%%%%
%%%%%%%%%Բ����%%%%%%%%%%
da = acosd((2*R^2-d^2)/(2*R^2));
faimax = 360-Nx*da;
faimin = 0;
Dx = da*ones(Nx,Ny,NP);
Dx(1,:,:) = 0;
Dx = cumsum(Dx,1);
fx = rand(Nx,Ny,NP)*(faimax-faimin)+faimin;   
fx = sort(fx,1);
fx1 = Dx+fx;
%%%%%%%%%������%%%%%%%%%%
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
%%%%%%%%%��Ⱥ%%%%%%%%%%%
f = fy+i.*fx;
ff = fy1+i.*fx1;
for g = 1:G
    g
   %%%%%%������Ӧ�ȣ�����ֵ�԰�%%%%
    for i = 1:NP
        Fit(i) = func_sparse_cylinder(k,NE,NA,seta0,fai0,ff(:,:,i));
    end
    maxFit = max(Fit);           %���ֵ
    minFit = min(Fit);           %��Сֵ
    rr = find(Fit==maxFit);      %�ҳ����ֵ
    fBest = ff(:,:,rr(1,1));     %�������Ÿ���
    Fit = (Fit-minFit)/(maxFit-minFit);  %��һ����Ӧ��ֵ
    %%%%%%%�������̶ĵĸ��Ʋ���%%%%%%
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
    %%%%%%%%%%���ڸ��ʵĽ������%%%%%%%%
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
    %%%%%%%���ڸ��ʵı������%%%%%%%
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
    ff(:,:,1) = fBest;                %�������Ÿ���������Ⱥ��
    trace(g) = maxFit;                %����������Ӧ��
end
figure
plot(trace)
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')
grid on
save fBest.mat fBest