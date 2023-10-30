%%%%%%%%%%%%%%%%¼ÆËã×î´óÅÔ°ê%%%%%%%%%%%%%%%
function   MSLL = func_cylinder(k,NE,NA,M,N,seta0,fai0,f0)
bottom = -50;
f = reshape(f0,N,M);
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
    F1 = exp(-sqrt(-1)*2*pi*rr.*(sin(seta2).*(sin(seta0)...
        *cos(fai0-faimn)-sin(seta(ne))*cos(fai0-faimn))...
        +cos(seta2)*(cos(seta0)-cos(seta(ne))))).*f;
    F(ne) = abs(sum(sum(F1)));
end
FdB1 = 20*log10(F/max(max(F)));
mm = ceil(find(FdB1==max(FdB1)));
tu_up = 0;
while (FdB1(mm+tu_up)>=FdB1(mm+tu_up+1))
    tu_up = tu_up+1;
end
tu_down = 0;
while (FdB1(mm-tu_down)>=FdB1(mm-tu_down-1))
    tu_down = tu_down+1;
end
FdB1(mm-tu_down:mm+tu_up) = bottom;
sll_1 = max(FdB1);
for na = 1:NA
    F1 = exp(-sqrt(-1)*2*pi*rr.*(sin(seta2).*(sin(seta0)...
        *cos(fai0-faimn)-sin(seta0)*cos(fai(na)-faimn))...
        +cos(seta2)*(cos(seta0)-cos(seta0)))).*f;
    F(na) = abs(sum(sum(F1)));
end
FdB2 = 20*log10(F/max(max(F)));
nn = find(FdB2==max(FdB2));
tv_up = 0;
while (FdB2(nn+tv_up)>=FdB2(nn+tv_up+1))
    tv_up = tv_up+1;
end
tv_down = 0;
while (FdB2(nn-tv_down)>=FdB2(nn-tv_down-1))
    tv_down = tv_down+1;
end
FdB2(nn-tv_down:nn+tv_up) = bottom;
sll_2 = max(FdB2);
MSLL = min(abs(sll_1),abs(sll_2));