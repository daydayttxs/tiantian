function MSLL = func_line(d,lamda,NN,seta0,f)
M = length(f);
bottom = -50;
seta = linspace(-pi/2,pi/2,NN);
for m=1:NN
    fai = 2*pi*d/lamda*(0:(M-1))*(sin(seta(m))-seta0);
    F1(m) = abs(sum(exp(sqrt(-1)*fai).*f'));
end
%%%%%%%%%%%% Çó×î´óÅÔ°ê  %%%%%%%%
FdB = 20*log10(F1/max(F1));
maxFdB = max(FdB);
rr=find(FdB == maxFdB);
mm = rr(1);
tu_up = 0;
while (FdB(mm+tu_up)>=FdB(mm+tu_up+1))
    tu_up = tu_up+1;
end
tu_down = 0;
while (FdB(mm-tu_down)>=FdB(mm-tu_down-1))
    tu_down = tu_down+1;
end
FdB(mm-tu_down:mm+tu_up) = bottom;
MSLL = max(FdB);
