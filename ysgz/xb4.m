%  1-D�ź�ѹ�����е�ʵ��(����ƥ��׷�ٷ�Orthogonal Matching Pursuit)
%  ������M>=K*log(N/K),K��ϡ���,N�źų���,���Խ�����ȫ�ع�
%  �����--��۴�ѧ���ӹ���ϵ ɳ��  Email: wsha@eee.hku.hk

clc;clear
%%  1. ʱ������ź�����
K=100;      %  ϡ���(��FFT���Կ�����)
N=2048;    %  �źų���
M=614;     %  ������(M>=K*log(N/K),����40,���г���ĸ���)

load('image_data2point.mat')
A=image_data;
B=A.';
x=B(64,1:2048);
%%  2.  ʱ���ź�ѹ������
Phi=randn(M,N);                                   %  ��������(��˹�ֲ�������)
s=Phi*x.';                                        %  ������Բ��� 

%%  3.  ����ƥ��׷�ٷ��ع��ź�(��������L_1�������Ż�����)
m=2*K;    %  �㷨��������(m>=K)
wtype = 'sym5';
wlev_max = wmaxlev(N,wtype);
if wlev_max == 0
    fprintf('\nThe parameter N and wtype does not match!\n');
end
dwtmode('per');
for wlev = 1:wlev_max
    ww = dwtmtx(N,wtype,wlev);
%     x = randn(1,N);
%     y1 = (ww*x')';
% Psi=ww;
% Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����                                   
end
Psi=ww;
T=Phi*Psi';                                   %  �ָ�����(��������*�������任����)
hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=s;                                            %  �в�ֵ

for times=1:m;                                    %  ��������(�������������,�õ�������ΪK)
    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������
    T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  ��С����,ʹ�в���С
    r_n=s-Aug_t*aug_y;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
end
hat_y(pos_array)=aug_y;                           %  �ع�����������
hat_x=real(Psi'*hat_y.');                         %  ���渵��Ҷ�任�ع��õ�ʱ���ź�

%%  4.  �ָ��źź�ԭʼ�źŶԱ�
figure(1);
hold on;
plot(hat_x,'k.-')                                 %  �ؽ��ź�
plot(x,'r')                                       %  ԭʼ�ź�
legend('�ָ��ź�','ԭʼ�ź�')
f=norm(hat_x.'-x)/norm(x)                           %  �ع����
a=sum((x'-hat_x).^2)/N;
err=sqrt(a)/max(x');
