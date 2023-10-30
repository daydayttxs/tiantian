%%%%%%%%%%%%%%%%¼ÆËã×î´óÅÔ°ê%%%%%%%%%%%%%%%
function   MSLL = func_sparse_cylinder(k,NE,NA,seta0,fai0,f0)
bottom = -50;
h = real(f0)+eps;
rr = sqrt(k^2+h.^2);
faimn = imag(f0)/180*pi;
seta1 = atan(k./h);                             
fai = linspace(0,2*pi,NA);                                       
seta = linspace(0,pi,NE); 
for ne = 1:NE
    F1 = exp(-sqrt(-1)*2*pi*rr.*(sin(seta1).*(sin(seta0)*...
        cos(fai0-faimn)-sin(seta(ne))*cos(fai0-faimn))...
        +cos(seta1)*(cos(seta0)-cos(seta(ne)))));
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
    F1 = exp(-sqrt(-1)*2*pi*rr.*(sin(seta1).*(sin(seta0)*...
        cos(fai0-faimn)-sin(seta0)*cos(fai(na)-faimn))...
        +cos(seta1)*(cos(seta0)-cos(seta0))));
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