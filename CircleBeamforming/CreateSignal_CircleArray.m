function [signal,len_num,num]=CreateSignal_CircleArray(style,signal_scan,R,N,f0,fs,T,alfa,B,input)
%--------------2012-02-10------------------------------
% ����Բ���������ź�
% ����:        style ���ź���ʽ     1---խ���źţ���Ƶ�źţ�    2---����źţ�LFM��  3---ָ���ź�
%        signal_scan : �ɽ����źŵĿ��Ƿ�Χ ����λ���ȣ�
%                  R : Բ��뾶
%                  N : ��Ԫ����
%                 f0 : խ����Ƶ�ʻ��߿��������Ƶ��
%                 fs ��������
%                  T ���ź�ʱ��
%               alfa : ���䷽�򣨵�λ���ȣ�
%                  B ���źŴ���
%              input : ��style=3��ʱ��ò�����Ч�������ָ���ź�
% �����      signal ��Բ������븴�ź�
%            len_num : �������źŵ���Ԫ��Ŀ.
%                num : �������źŵ���Ԫ���
%
% �㷨�� ����խ���ĵ�Ƶ�źŲ������Ƶķ����γ������ź�
%        ���ڿ��LFM�źŲ���ʱ�ӷ����γ������ź�
%---------------------------------------------------------
% ��ע
% 2014-02-28�޸ģ�����style=3�����������֪�źŵ������ͨ��Ϊͨ���źţ�
c=1500;%%Ĭ������Ϊ1500m/s
thita=alfa*pi/180;
signal_scan=signal_scan*pi/180;
signal=0;
%%--ȷ�����źŵ���Ԫ---
scan_begin = thita-signal_scan/2;
if(scan_begin<=0)
    scan_begin = mod(scan_begin + 2*pi,2*pi);
end
scan_end = thita + signal_scan/2;
if(scan_end > 2*pi)
    scan_end = mod( scan_end-2*pi,2*pi);
end
kk=1;
for k = 1:N   %���ڼ������Ԫ�ǣ�Ĭ��180�ȵĿ��Ƿ�Χ�ڵ���Ԫ�����Խӵ��ź�
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
len_num=length(num); %%�����������Ԫ���Խ��յ��źţ����ź���Բ�ĺ���Ԫ���ߵķ��������䣬��֮��ż������Ԫ���Խӵ��źš�
% polar(psi,p1,'mo'), hold on; %ͼʾ��Ԫ
if(style==1)                    %%%-------- խ���ź����� -----------
    tao0=zeros(1,N);
    sig=zeros(N,T*fs);
    s=exp(j*2*pi*f0*(0:1/fs:T-1/fs));
    sig(num,:)=ones(len_num,1)*s;
    tao0(num)=exp(-j*2*pi*f0*R*(1-cos(thita-psi(num)))/c);
    tao=tao0.'*ones(1,T*fs);
    signal=sig.*tao;
elseif(style==2)                 %%%-------- ����ź����� --------------
    K=B/T;
    t=0:1/fs:T-1/fs;
    s=exp(j*(2*pi*(f0-B/2)*t+pi*K*t.^2));
    beishu=1;
    fs_new=beishu*fs;    
    s_new=resample(s,beishu,1);  %%�����ز���
    lens_new=length(s_new);
    MaxTao=round(R/c*fs_new); %%�����ʱ����ΪR��MaxTao�����ʱ����
    sig_ru=zeros(N,2*MaxTao+lens_new);  %%�ز�����������ź������ʼ��    
    Ntao=zeros(1,N);
    Ntao(num)=fix(R*(1-cos(thita-psi(num)))/c*fs_new);
    for k=1:len_num
        sig_ru(num(k),MaxTao+Ntao(num(k))+1:MaxTao+Ntao(num(k))+lens_new)=s_new;
    end
    sig_ru=sig_ru.';
    sig_ru=resample(sig_ru,1,beishu);
    signal=sig_ru.';  %%ÿһ�о���һ����Ԫ���յ��ź�
    sig_len=size(signal,2);
elseif(style==3)
    s=input;
    beishu=1;
    fs_new=beishu*fs;    
    s_new=resample(s,beishu,1);  %%�����ز���
    clear s;
    lens_new=length(s_new);
    MaxTao=round(R/c*fs_new); %%�����ʱ����ΪR��MaxTao�����ʱ����
    sig_ru=zeros(N,10+MaxTao+lens_new);  %%�ز�����������ź������ʼ��    
    Ntao=zeros(1,N);
    Ntao(num)=fix(R*(1-cos(thita-psi(num)))/c*fs_new);
    for k=1:len_num
        sig_ru(num(k),MaxTao+Ntao(num(k))+1:MaxTao+Ntao(num(k))+lens_new)=s_new;
    end
    sig_ru=sig_ru.';
    sig_ru=resample(sig_ru,1,beishu);
    signal=sig_ru.';  %%ÿһ�о���һ����Ԫ���յ��ź�
    sig_len=size(signal,2);
    
else
    error('�ź���������');
end




