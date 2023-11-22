%% 仿真：圆阵波束形成
%% 算法(method)：“1”---时域的时延波束形成
%%               “2”---相移波束形成
%%               “3”---频域波束形成
%% 信号(style)： “1”---窄带信号（单频信号）
%%               “2”---宽带信号 (LFM信号）
%%----------------------------------------------------------------------
clc
clear all
close all
style=1;
method=[3];
if (style==1) disp('信号类型：圆阵窄带信号');
else disp('信号类型：圆阵宽带信号');
end
lenm=length(method);
for k=1:lenm
    if (method(k)==1)
        disp('算法：时延波束形成')
    elseif(method(k)==2)
        disp('算法：相移波束形成')
    elseif(method(k)==3)
        disp('算法：频域波束形成');
    end
end
thita=90;  %%信号入射方位
signal_scan=360; %%生成信号时可接收信号的开角范围
N=64;
dthita=5.625;%%波束形成时，搜索角度旋转步长
scan=120;%%波束形成时的开角范围
disp(['仿真信号入射方位：' num2str(thita) '度']);
disp(['制造信号的开角范围：' num2str(signal_scan) '度']);
disp(['搜索角度旋转步长：' num2str(dthita) '度']);
disp(['波束形成时的开角范围：' num2str(scan) '度']);

c=1500;
f0=3000;
fs=24000%%24000;%20*f0;
T=0.05%%0.34;%;
B=2000;
R=1.291;%(N-2)*c/(4*pi*f0);%1.1937%%
[signal,num]=CreateSignal_CircleArray(style,signal_scan,R,N,f0,fs,T,thita,B);
%%------ 看每个阵元上的信号波形-----------
figure(1);
for k=1:N
    plot(real(signal(k,1:length(signal)/10))+(k-1)*2,'*-r');hold on;
end
grid on;title('每个阵元上的信号波形');
%%-----------------------------------------
[power,powerdB,doa,sigout1]=Beamforming_CircleArray(style,method,scan,dthita,signal,R,N,f0,fs);
figure(2);plot([0:dthita:359],power,'b');grid on;hold on;xlabel('方位/度');ylabel('能量');
figure(3);plot([0:dthita:359],powerdB,'*-b');grid on;hold on;
xlabel('方位/度');ylabel('能量/分贝');
Norm=0.5*max(abs(powerdB));
for k=1:length(powerdB)
    if powerdB(k)<-Norm
        powerdB(k)=-Norm;
    end
     powerdB(k)= powerdB(k)+Norm;
end
figure(4);polar([0:dthita:359]*pi/180,powerdB), hold on; %极坐标
disp(['波达方位估计为：' num2str(doa) '度']);