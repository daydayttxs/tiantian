%%%%%%%%%%%%优化后方向图%%%%%%%%%%%%%
clc;close all;clear all
load('fBest.mat')
f0 = fBest;
k = 2;
seta0 = pi/2;
fai0 = pi;
NE = 360;
NA = 360;
%%%%%%%%%%%%%%%%%%%%%%
h = real(f0)+eps;
rr = sqrt(k^2+h.^2);
faimn = imag(f0)/180*pi;
seta1 = atan(k./h);
fai = linspace(0,2*pi,NA);    
seta = linspace(0,pi,NE); 
for ne = 1:NE
    for na = 1:NA
        F1 = exp(-sqrt(-1)*2*pi*rr.*(sin(seta1).*(sin(seta0)...
            *cos(fai0-faimn)-sin(seta(ne))*cos(fai(na)-faimn))...
            +cos(seta1)*(cos(seta0)-cos(seta(ne)))));
        F(na,ne) = abs(sum(sum(F1)));
    end
end
figure           %三维立体图
[X,Y] = meshgrid((fai*180/pi),(seta*180/pi));
Z = 20*log10(F/max(max(F)));
number = find(Z<-50);
g_temp = -50+unifrnd(-1,1,1,length(number));
for ii = 1:length(number);
    Z(number(ii)) = g_temp(ii);
end
mesh(X,Y,Z)
xlabel('\theta俯仰角（度）')
ylabel('\phi方位角（度）')
zlabel('阵列增益(dB)')
figure %方位向切面图
temp1 = Z(:,round(NE*((pi-seta0)/pi)));
plot(fai*180/pi,temp1)
grid
xlabel('\phi方位角（度）')
ylabel('阵列增益(dB)')
figure %俯仰向切面图
temp2 = Z(round(NA*((2*pi-fai0)/(2*pi))),:);
plot(seta*180/pi,temp2)
grid
xlabel('\theta俯仰角（度）')
ylabel('阵列增益(dB)')
figure
X=real(fBest);
Y=imag(fBest);
plot(Y,X,'*b')
xlabel('圆环向')
ylabel('俯仰向')