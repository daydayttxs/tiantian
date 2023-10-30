%%%%%%%%%%����Բ����%%%%%%%%%%
%%%%%%%%%%%%��ʼ������%%%%%%%%%%%%
clear all;                %������б���
close all;                %��ͼ
clc;                      %����
M = 20;                   %����������Բ����
N = 24;                   %ÿ��Բ������Ԫ��
k = 2;                    %��Ԫ�뾶�벨��֮��      
seta0 = pi/2;             %�źŵ�����λ��           
fai0 = pi;                %�źŵ���������
NA = 360;                 %�ռ䷽λ�ǲ�����
NE = 360;                 %�ռ丩���ǲ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 0:0.5:0.5*(M-1);      
r = sqrt(k^2+h.^2) ;
seta1 = zeros(1,M);
seta1(1) = pi/2;
seta1(2:M) = atan(k./h(2:M));
fai = linspace(0,2*pi,NA);   
seta = linspace(0,pi,NE); 
fain = (0:(N-1))*2*pi/N;
faimn = repmat(fain',1,M);
rr = repmat(r,N,1);
seta2 = repmat(seta1,N,1);
for ne = 1:NE
    for na = 1:NA
        F1 = exp(-sqrt(-1)*2*pi*rr.*(sin(seta2).*(sin(seta0)...
            *cos(fai0-faimn)-sin(seta(ne))*cos(fai(na)-faimn))...
            +cos(seta2)*(cos(seta0)-cos(seta(ne)))));
        F(na,ne) = abs(sum(sum(F1)));
    end
end
figure           %��ά����ͼ
[X,Y] = meshgrid((seta*180/pi),(fai*180/pi));
Z = 20*log10(F/max(max(F)));
number = find(Z<-50);
g_temp = -50+unifrnd(-1,1,1,length(number));
for ii = 1:length(number);
    Z(number(ii)) = g_temp(ii);
end
mesh(X,Y,Z)
xlabel('\theta�����ǣ��ȣ�')
ylabel('\phi��λ�ǣ��ȣ�')
zlabel('��������(dB)')
figure %��λ������ͼ
temp1 = Z(:,round(NE*((pi-seta0)/pi)));
plot(fai*180/pi,temp1)
grid
xlabel('\phi��λ�ǣ��ȣ�')
ylabel('��������(dB)')
figure %����������ͼ
temp2 = Z(round(NA*((2*pi-fai0)/(2*pi))),:);
plot(seta*180/pi,temp2)
grid
xlabel('\theta�����ǣ��ȣ�')
ylabel('��������(dB)')