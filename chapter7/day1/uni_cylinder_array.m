%%%%%%%%%%均匀圆柱阵%%%%%%%%%%
%%%%%%%%%%%%初始化参数%%%%%%%%%%%%
clear all;                %清除所有变量
close all;                %清图
clc;                      %清屏
M = 3;                   %拄面阵列中圆环数
N = 256;                   %每个圆环中阵元数
k = 4.8;                    %阵元半径与波长之比      
seta0 = pi/2;             %信号到来方位角           
fai0 = pi;                %信号到来俯仰角
NA = 360;                 %空间方位角采样数
NE = 360;                 %空间俯仰角采样数
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
figure           %三维立体图
[X,Y] = meshgrid((seta*180/pi),(fai*180/pi));
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