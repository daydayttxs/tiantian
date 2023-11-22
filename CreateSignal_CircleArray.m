function [signal,len_num,num]=CreateSignal_CircleArray(style,signal_scan,R,N,f0,fs,T,alfa,B,input)
%--------------2012-02-10------------------------------
% 制造圆柱阵入射信号
% 输入:        style ：信号形式     1---窄带信号（单频信号）    2---宽带信号（LFM）  3---指定信号
%        signal_scan : 可接收信号的开角范围 （单位：度）
%                  R : 圆阵半径
%                  N : 阵元个数
%                 f0 : 窄带的频率或者宽带的中心频率
%                 fs ：采样率
%                  T ：信号时长
%               alfa : 入射方向（单位：度）
%                  B ：信号带宽
%              input : 当style=3的时候该参量有效，输入的指定信号
% 输出：      signal ：圆阵的输入复信号
%            len_num : 有输入信号的阵元数目.
%                num : 有输入信号的阵元序号
%
% 算法： 对于窄带的单频信号采用相移的方法形成输入信号
%        对于宽带LFM信号采用时延方法形成输入信号
%---------------------------------------------------------
% 备注
% 2014-02-28修改：加入style=3的情况：给已知信号的情况（通常为通信信号）
c=1500;%%默认声速为1500m/s
thita=alfa*pi/180;
signal_scan=signal_scan*pi/180;
signal=0;
%%--确定有信号的阵元---
scan_begin = thita-signal_scan/2;
if(scan_begin<=0)
    scan_begin = mod(scan_begin + 2*pi,2*pi);
end
scan_end = thita + signal_scan/2;
if(scan_end > 2*pi)
    scan_end = mod( scan_end-2*pi,2*pi);
end
kk=1;
for k = 1:N   %用于计算的阵元角，默认180度的开角范围内的阵元都可以接到信号
    p1(k) = 0;
    psi(k) = 2 * pi *(k-1) / N;
    if( scan_begin >= scan_end)
        if(scan_begin / (psi(k)+eps) <= 1.0+eps)
            p1(k) = R ;  num(kk)=k;  kk=kk+1;
        elseif( psi(k) / (scan_end+eps) <= 1.0+eps)
            p1(k) = R ;  num(kk)=k;  kk=kk+1;
        end
    else
        if ( scan_begin / (psi(k)+eps) <= 1.0+eps) &  (psi(k) / (scan_end+eps) <= 1.0+eps)
            p1(k) = R;  num(kk)=k;  kk=kk+1;
        end
    end
end
len_num=length(num); %%如果奇数个阵元可以接收到信号，即信号在圆心和阵元连线的方向上入射，反之，偶数个阵元可以接到信号。
% polar(psi,p1,'mo'), hold on; %图示阵元
if(style==1)                    %%%-------- 窄带信号入射 -----------
    tao0=zeros(1,N);
    sig=zeros(N,T*fs);
    s=exp(j*2*pi*f0*(0:1/fs:T-1/fs));
    sig(num,:)=ones(len_num,1)*s;
    tao0(num)=exp(-j*2*pi*f0*R*(1-cos(thita-psi(num)))/c);
    tao=tao0.'*ones(1,T*fs);
    signal=sig.*tao;
elseif(style==2)                 %%%-------- 宽带信号入射 --------------
    K=B/T;
    t=0:1/fs:T-1/fs;
    s=exp(j*(2*pi*(f0-B/2)*t+pi*K*t.^2));
    beishu=1;
    fs_new=beishu*fs;    
    s_new=resample(s,beishu,1);  %%按列重采样
    lens_new=length(s_new);
    MaxTao=round(R/c*fs_new); %%最大延时长度为R，MaxTao最大延时点数
    sig_ru=zeros(N,2*MaxTao+lens_new);  %%重采样后的入射信号数组初始化    
    Ntao=zeros(1,N);
    Ntao(num)=fix(R*(1-cos(thita-psi(num)))/c*fs_new);
    for k=1:len_num
        sig_ru(num(k),MaxTao+Ntao(num(k))+1:MaxTao+Ntao(num(k))+lens_new)=s_new;
    end
    sig_ru=sig_ru.';
    sig_ru=resample(sig_ru,1,beishu);
    signal=sig_ru.';  %%每一行就是一个阵元接收的信号
    sig_len=size(signal,2);
elseif(style==3)
    s=input;
    beishu=1;
    fs_new=beishu*fs;    
    s_new=resample(s,beishu,1);  %%按列重采样
    clear s;
    lens_new=length(s_new);
    MaxTao=round(R/c*fs_new); %%最大延时长度为R，MaxTao最大延时点数
    sig_ru=zeros(N,10+MaxTao+lens_new);  %%重采样后的入射信号数组初始化    
    Ntao=zeros(1,N);
    Ntao(num)=fix(R*(1-cos(thita-psi(num)))/c*fs_new);
    for k=1:len_num
        sig_ru(num(k),MaxTao+Ntao(num(k))+1:MaxTao+Ntao(num(k))+lens_new)=s_new;
    end
    sig_ru=sig_ru.';
    sig_ru=resample(sig_ru,1,beishu);
    signal=sig_ru.';  %%每一行就是一个阵元接收的信号
    sig_len=size(signal,2);
    
else
    error('信号类型有误');
end




