function [power,powerdB,doa,signal_out]=Beamforming_CircleArray(style,method,scan,dthita,sig,R,N,f0,fs,W)
%--------------2012-02-17------------------------------
% 圆弧阵波束形成，360度全方位扫描
% 输入:        style ：信号形式     1---窄带信号（单频信号）    2---宽带信号（LFM）
%             method : 波束形成算法类型  
%                                 “1”--- 时域的时延波束形成
%                                 “2”--- 相移波束形成
%                                 “3”--- 频域波束形成
%                                 “4”--- 经过三次FFT的频域波束形成，延迟滤波器和输入信号在频域相乘
%              scan : 开角范围（单位：度）
%           dthita : 扫描步长 （单位：度）
%                sig : 圆阵输入信号
%                  R : 圆阵半径
%                  N : 阵元个数
%                 f0 : 窄带的中心频率或者宽带的最小频率
%                 fs ：采样率
%                  W : 加权值
% 输出：       power ：能量
%            powerdB : 能量归一化分贝输出
%              doa : 信号波达方向( 单位：度）
%         signal_out : 波束在预成方向上叠加后的输出信号，每一行代表一个方向上的输出
%---------------------------------------------------------
c=1500;  %%默认声速为1500m/s
[sig_N,sig_len]=size(sig);
if( sig_N ~= N )
    error('信号和圆阵个数不匹配');
end
if ( style ~=1 & style ~=2  )
    error('请正确选择输入信号类型');
end
if ( method ~=1 & method ~=2 & method ~=3 & method ~=4 )
    error('请正确选择波束形成算法');
end

psi = 2 * pi *(0:N-1) / N; %%圆心角
thita=[0:dthita:359]*pi/180;
power=zeros(1,length(thita));
scan_begin = thita - scan *pi/180/ 2;
scan_begin( scan_begin < 0 ) = scan_begin( scan_begin < 0 )+2*pi;
scan_end = thita + scan * pi/180/ 2;
scan_end(scan_end >= 2*pi) = scan_end(scan_end >= 2*pi)-2*pi;
flag=zeros(length(thita),N); %% 标志位，对应位为1时，表示该号阵元参与计算
%%----------  确定预成方向上参与计算的阵元号,存入flag中  ----
for a=thita
            index=round(a/(dthita*pi/180)+1);
            if( scan_end (index) / (scan_begin (index)+eps) <= (1.0+0.0001))
                flag(index,find(scan_begin (index)./ (psi+eps) <=(1.0+eps)))=1;
                flag(index,find((psi-eps) ./(scan_end (index)+eps) <= (1.0+eps)))=1; 
            else
                flag(index,find(scan_begin (index)./(psi+eps) <= (1.0+eps) & psi ./(scan_end (index)+eps)<= (1.0+eps)))=1;
            end 
end
%%-------------------------------------------------------------
if( method == 1)                %%----  时延波束形成对入射信号不限制，窄带、宽带均可  ( method == 1) ---- 
    beishu=1;
    %%---- 升采样（可选） ----
%         beishu=10;
%         sig=resample(sig.',beishu,1);
%         sig=sig.';
%         fs=fs*beishu;
%         sig_len=sig_len*beishu;
    %%------------------------
    MaxTao=round(2*R/c*fs);
    for a=thita
%         a*180/pi
        sig_buchang=zeros(N,sig_len+MaxTao);
        index=round(a/(dthita*pi/180)+1);
        Number_Array_InUse=find(flag(index,:)==1);
        Ntao0=round(R*(1-cos(a-psi(Number_Array_InUse)))/c*fs);
        for kk=1:length(Number_Array_InUse)
            sig_buchang(Number_Array_InUse(kk),MaxTao-Ntao0(kk)+1:MaxTao-Ntao0(kk)+sig_len)=sig(Number_Array_InUse(kk),:);
        end
        sig_buchang=resample(sig_buchang.',1,beishu);
        sig_buchang=sig_buchang.';
        power(index)=sum(sum(real(sig_buchang(Number_Array_InUse,:))).^2);
        signal_out(index,:)=sum(sig_buchang(Number_Array_InUse,:));   %%注意：这里的输出的波束信号有冗余
    end

elseif( method == 2)            %%----  相移波束形成只针对窄带  ( method == 2)  ----  
    if( style~=1 )
        error('相移波束形成只针对窄带，请正确输入信号类型');
    else        
        for a=thita 
            index=round(a/(dthita*pi/180)+1);
            Number_Array_InUse=find(flag(index,:)==1);
            tao0=exp(1i*2*pi*f0*R*(1-cos(a-psi(Number_Array_InUse)))/c);
            tao_array=tao0.'*ones(1,sig_len);
            sig_buchang=sig(Number_Array_InUse,:).*tao_array;
            power(index)=sum(sum(real(sig_buchang)).^2);
            signal_out(index,:)=sum(sig_buchang);
        end   
    end

elseif( method == 3 )                       %%----  频域波束形成对入射信号不限制，窄带、宽带均可   ( method == 3) ----
    sig=sig.';
    fft_sig=fft(sig);
    fft_sig=fft_sig.';
    fft_index=(0:sig_len-1)/sig_len*fs;
    for a=thita
        index=round(a/(dthita*pi/180)+1);
        Number_Array_InUse=find(flag(index,:)==1);
        tao0=-R*(1-cos(a-psi(Number_Array_InUse)))/c;
        tao=tao0.'*ones(1,floor(sig_len/2)-2);  %%length(Number_Array_InUse)行，每行中的元素相同，该阵元号码的延时
        freq=repmat(fft_index(2:floor(sig_len/2)-1),length(Number_Array_InUse),1);
        A=exp(-1i*2*pi*freq.*tao);
        freq_buchang=2*fft_sig(Number_Array_InUse,2:floor(sig_len/2)-1).*A;
        ifft_fre=ifft(sum(freq_buchang));
        power(index)=sum(real(ifft_fre).^2);   
        signal_out(index,:)=ifft_fre;
    end
elseif( method == 4 )  %%----  三次FFT的频域波束形成对入射信号不限制，窄带、宽带均可   ( method == 4) ----
    if(nargin==9) %%没有加权值时
        W = ones(sig_len,length(thita));
    end
    %%-------1. 将各个阵元接收到的信号变换到频域---------------
    fft_sig=fft(sig,[],2);                        %%第一次FFT
    disp('第一次FFT');
    %%-------2.将输入信号各个频点上对应的值进行FFT变换---------------
    K=100;                                         %%一次处理K点
    NFFT=2^nextpow2(K+N-1);                        %%每次傅里叶变换的点数
    fft_sig=fft_sig(end:-1:1,:);                   %%注意，这里要取反一次才能保证结果正确
    fft_sig=[fft_sig;zeros(NFFT-N,sig_len)];       %%将fft_sig的每一列都补成NFFT点 
    disp('第二次FFT');
    fft2_sig=fft(fft_sig,NFFT);                    %%fft2_sig的每一列代表一个频点上所有N点的fft，行数代表共有多少频点  %%第二次FFT
    %%-------3.将延迟滤波器传输函数变换到频域---------------
    alfa=repmat(thita',1,N)-repmat(psi,length(thita),1);             %%矩阵大小length(thita)*N
    alfa=alfa.*flag;                               %%将接收开角信息加入，但是还没有利用上阵列流型的对称性
    alfa=2*pi*R*(1-cos(alfa))/c;
    alfa=reshape(alfa',1,length(thita)*N);    
    fk=(0:sig_len-1)/sig_len*fs;
    for k1=1:sig_len
        beta=exp(1i*fk(k1)*alfa);
        beta=[beta zeros(1,K-mod(length(thita)*N,K))];
        fa=[];
        for k2=1:K:length(thita)*N
            if k2==1
                fa((k2-1)/K+1,:)=[zeros(1,N-1),beta(k2:k2+K-1)];        %%第一次补零点
            else
                fa((k2-1)/K+1,:)=[beta(k2-(N-1):k2-1) beta(k2:k2+K-1)];  %%在每一段的前边补上前一段保留下来的(M-1)点               
            end
        end
        fa=[fa zeros(size(fa,1),NFFT-size(fa,2))];   %%将a的每一行都补成NFFT点，a的行数代表截断的次数，a的列数代表一次处理的阵列流型的数据
        fft_a=fft(fa,NFFT,2);                        %%第三次FFT
        fk_a=fft_a;
        sigfk=fft2_sig(:,k1);
        s=repmat(sigfk.',size(fft_a,1),1);
        Y=s.*fk_a;
        y=ifft(Y,NFFT,2);
        r=y(:,N-1+1:N+K-1);
        out=reshape(r.',1,size(r,1)*size(r,2));
        rr(k1,:)=out(N:N:length(thita)*N);
    end
    rr=rr.*W;
    power=sum(abs(rr));
    signal_out=(ifft(rr)).';
end
powerdB=10*log10((power+eps)/(max(abs(power))+eps));
doa=(find(powerdB==max(powerdB))-1)*dthita;



        
        
        
    