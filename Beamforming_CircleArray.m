function [power,powerdB,doa,signal_out]=Beamforming_CircleArray(style,method,scan,dthita,sig,R,N,f0,fs,W)
%--------------2012-02-17------------------------------
% Բ�������γɣ�360��ȫ��λɨ��
% ����:        style ���ź���ʽ     1---խ���źţ���Ƶ�źţ�    2---����źţ�LFM��
%             method : �����γ��㷨����  
%                                 ��1��--- ʱ���ʱ�Ӳ����γ�
%                                 ��2��--- ���Ʋ����γ�
%                                 ��3��--- Ƶ�����γ�
%                                 ��4��--- ��������FFT��Ƶ�����γɣ��ӳ��˲����������ź���Ƶ�����
%              scan : ���Ƿ�Χ����λ���ȣ�
%           dthita : ɨ�貽�� ����λ���ȣ�
%                sig : Բ�������ź�
%                  R : Բ��뾶
%                  N : ��Ԫ����
%                 f0 : խ��������Ƶ�ʻ��߿������СƵ��
%                 fs ��������
%                  W : ��Ȩֵ
% �����       power ������
%            powerdB : ������һ���ֱ����
%              doa : �źŲ��﷽��( ��λ���ȣ�
%         signal_out : ������Ԥ�ɷ����ϵ��Ӻ������źţ�ÿһ�д���һ�������ϵ����
%---------------------------------------------------------
c=1500;  %%Ĭ������Ϊ1500m/s
[sig_N,sig_len]=size(sig);
if( sig_N ~= N )
    error('�źź�Բ�������ƥ��');
end
if ( style ~=1 & style ~=2  )
    error('����ȷѡ�������ź�����');
end
if ( method ~=1 & method ~=2 & method ~=3 & method ~=4 )
    error('����ȷѡ�����γ��㷨');
end

psi = 2 * pi *(0:N-1) / N; %%Բ�Ľ�
thita=[0:dthita:359]*pi/180;
power=zeros(1,length(thita));
scan_begin = thita - scan *pi/180/ 2;
scan_begin( scan_begin < 0 ) = scan_begin( scan_begin < 0 )+2*pi;
scan_end = thita + scan * pi/180/ 2;
scan_end(scan_end >= 2*pi) = scan_end(scan_end >= 2*pi)-2*pi;
flag=zeros(length(thita),N); %% ��־λ����ӦλΪ1ʱ����ʾ�ú���Ԫ�������
%%----------  ȷ��Ԥ�ɷ����ϲ���������Ԫ��,����flag��  ----
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
if( method == 1)                %%----  ʱ�Ӳ����γɶ������źŲ����ƣ�խ�����������  ( method == 1) ---- 
    beishu=1;
    %%---- ����������ѡ�� ----
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
        signal_out(index,:)=sum(sig_buchang(Number_Array_InUse,:));   %%ע�⣺���������Ĳ����ź�������
    end

elseif( method == 2)            %%----  ���Ʋ����γ�ֻ���խ��  ( method == 2)  ----  
    if( style~=1 )
        error('���Ʋ����γ�ֻ���խ��������ȷ�����ź�����');
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

elseif( method == 3 )                       %%----  Ƶ�����γɶ������źŲ����ƣ�խ�����������   ( method == 3) ----
    sig=sig.';
    fft_sig=fft(sig);
    fft_sig=fft_sig.';
    fft_index=(0:sig_len-1)/sig_len*fs;
    for a=thita
        index=round(a/(dthita*pi/180)+1);
        Number_Array_InUse=find(flag(index,:)==1);
        tao0=-R*(1-cos(a-psi(Number_Array_InUse)))/c;
        tao=tao0.'*ones(1,floor(sig_len/2)-2);  %%length(Number_Array_InUse)�У�ÿ���е�Ԫ����ͬ������Ԫ�������ʱ
        freq=repmat(fft_index(2:floor(sig_len/2)-1),length(Number_Array_InUse),1);
        A=exp(-1i*2*pi*freq.*tao);
        freq_buchang=2*fft_sig(Number_Array_InUse,2:floor(sig_len/2)-1).*A;
        ifft_fre=ifft(sum(freq_buchang));
        power(index)=sum(real(ifft_fre).^2);   
        signal_out(index,:)=ifft_fre;
    end
elseif( method == 4 )  %%----  ����FFT��Ƶ�����γɶ������źŲ����ƣ�խ�����������   ( method == 4) ----
    if(nargin==9) %%û�м�Ȩֵʱ
        W = ones(sig_len,length(thita));
    end
    %%-------1. ��������Ԫ���յ����źű任��Ƶ��---------------
    fft_sig=fft(sig,[],2);                        %%��һ��FFT
    disp('��һ��FFT');
    %%-------2.�������źŸ���Ƶ���϶�Ӧ��ֵ����FFT�任---------------
    K=100;                                         %%һ�δ���K��
    NFFT=2^nextpow2(K+N-1);                        %%ÿ�θ���Ҷ�任�ĵ���
    fft_sig=fft_sig(end:-1:1,:);                   %%ע�⣬����Ҫȡ��һ�β��ܱ�֤�����ȷ
    fft_sig=[fft_sig;zeros(NFFT-N,sig_len)];       %%��fft_sig��ÿһ�ж�����NFFT�� 
    disp('�ڶ���FFT');
    fft2_sig=fft(fft_sig,NFFT);                    %%fft2_sig��ÿһ�д���һ��Ƶ��������N���fft�����������ж���Ƶ��  %%�ڶ���FFT
    %%-------3.���ӳ��˲������亯���任��Ƶ��---------------
    alfa=repmat(thita',1,N)-repmat(psi,length(thita),1);             %%�����Сlength(thita)*N
    alfa=alfa.*flag;                               %%�����տ�����Ϣ���룬���ǻ�û���������������͵ĶԳ���
    alfa=2*pi*R*(1-cos(alfa))/c;
    alfa=reshape(alfa',1,length(thita)*N);    
    fk=(0:sig_len-1)/sig_len*fs;
    for k1=1:sig_len
        beta=exp(1i*fk(k1)*alfa);
        beta=[beta zeros(1,K-mod(length(thita)*N,K))];
        fa=[];
        for k2=1:K:length(thita)*N
            if k2==1
                fa((k2-1)/K+1,:)=[zeros(1,N-1),beta(k2:k2+K-1)];        %%��һ�β����
            else
                fa((k2-1)/K+1,:)=[beta(k2-(N-1):k2-1) beta(k2:k2+K-1)];  %%��ÿһ�ε�ǰ�߲���ǰһ�α���������(M-1)��               
            end
        end
        fa=[fa zeros(size(fa,1),NFFT-size(fa,2))];   %%��a��ÿһ�ж�����NFFT�㣬a����������ضϵĴ�����a����������һ�δ�����������͵�����
        fft_a=fft(fa,NFFT,2);                        %%������FFT
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



        
        
        
    