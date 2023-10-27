%%%%%%%%%%%%%%%%%优化后方向图%%%%%%
clc;close all;clear all;
load('fBest3.mat')
c=1540;
f=8.5e6;

lamda = c/f;
d = 0.6*lamda;
seta0=0*pi/180;
L = 64;
NN=1800;
%%%%%%%%%%%%%%%%%%%%%%
seta=linspace(-pi/2,pi/2,NN);
for m = 1: NN
    fai=2* pi *d/lamda*(0: (L-1))*(sin(seta(m))-seta0);
    F1(m) = abs(sum(exp(sqrt(-1)*fai).*fBest'));
end
FdB = 20*log10(F1/max(F1));
figure
plot(seta*180/pi,FdB)
xlabel('\theta/(°)')
ylabel('阵列增益/dB')
grid on
axis([-90,90,-60,0])
figure
plot(fBest,'.')
xlabel('阵元位置')
ylabel('阵元标识')
grid on
axis([1,L,0.5,1.5])