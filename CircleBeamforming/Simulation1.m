%% ���棺Բ�����γ�
%% �㷨(method)����1��---ʱ���ʱ�Ӳ����γ�
%%               ��2��---���Ʋ����γ�
%%               ��3��---Ƶ�����γ�
%% �ź�(style)�� ��1��---խ���źţ���Ƶ�źţ�
%%               ��2��---����ź� (LFM�źţ�
%%----------------------------------------------------------------------
clc
clear all
close all
style=1;
method=[3];
if (style==1) disp('�ź����ͣ�Բ��խ���ź�');
else disp('�ź����ͣ�Բ�����ź�');
end
lenm=length(method);
for k=1:lenm
    if (method(k)==1)
        disp('�㷨��ʱ�Ӳ����γ�')
    elseif(method(k)==2)
        disp('�㷨�����Ʋ����γ�')
    elseif(method(k)==3)
        disp('�㷨��Ƶ�����γ�');
    end
end
thita=90;  %%�ź����䷽λ
signal_scan=360; %%�����ź�ʱ�ɽ����źŵĿ��Ƿ�Χ
N=64;
dthita=5.625;%%�����γ�ʱ�������Ƕ���ת����
scan=120;%%�����γ�ʱ�Ŀ��Ƿ�Χ
disp(['�����ź����䷽λ��' num2str(thita) '��']);
disp(['�����źŵĿ��Ƿ�Χ��' num2str(signal_scan) '��']);
disp(['�����Ƕ���ת������' num2str(dthita) '��']);
disp(['�����γ�ʱ�Ŀ��Ƿ�Χ��' num2str(scan) '��']);

c=1500;
f0=3000;
fs=24000%%24000;%20*f0;
T=0.05%%0.34;%;
B=2000;
R=1.291;%(N-2)*c/(4*pi*f0);%1.1937%%
[signal,num]=CreateSignal_CircleArray(style,signal_scan,R,N,f0,fs,T,thita,B);
%%------ ��ÿ����Ԫ�ϵ��źŲ���-----------
figure(1);
for k=1:N
    plot(real(signal(k,1:length(signal)/10))+(k-1)*2,'*-r');hold on;
end
grid on;title('ÿ����Ԫ�ϵ��źŲ���');
%%-----------------------------------------
[power,powerdB,doa,sigout1]=Beamforming_CircleArray(style,method,scan,dthita,signal,R,N,f0,fs);
figure(2);plot([0:dthita:359],power,'b');grid on;hold on;xlabel('��λ/��');ylabel('����');
figure(3);plot([0:dthita:359],powerdB,'*-b');grid on;hold on;
xlabel('��λ/��');ylabel('����/�ֱ�');
Norm=0.5*max(abs(powerdB));
for k=1:length(powerdB)
    if powerdB(k)<-Norm
        powerdB(k)=-Norm;
    end
     powerdB(k)= powerdB(k)+Norm;
end
figure(4);polar([0:dthita:359]*pi/180,powerdB), hold on; %������
disp(['���﷽λ����Ϊ��' num2str(doa) '��']);