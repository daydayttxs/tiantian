%%
clear ;clc
load C:\Users\Administrator\Documents\Vantage-4.3.0-2009141600\wy\data\d25mm_yuanxing\tiqu_data\RFData_SA1
% load E:\user-MATLAB\Test\实测\256合成孔径间隔1\data\提取数据\RFData_sig_water
a=256;   %%修改
%%
% % sig_dB=20*log10(abs(RFData_sig));
% % sig_water_dB=20*log10(abs(RFData_sig_water));
% % sig=sig_dB-sig_water_dB;
% sig=RFData_sig-RFData_sig_water;                                                                                                                                                                                                                                                                                                                                                                    
% sig1=sig;
% sig1(sig1==-Inf)=0;
% sig1(sig1==Inf)=0;
% sig1(sig1==0)=0;
% RFData_sig=sig1;
% % sig1=sig1-min(min(sig1));
RFData_sig=RFData_SA1;
RFData_sig(1:300,:)=RFData_sig(301:600,:);
% % RFData_sig_4(2500:3200,:)=0;
%%
env_dB=20*log10(abs(RFData_sig));
% env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
%% 信号截取
sum_a_data_need=env_dB';
sum_a_data_need_filter=sum_a_data_need;
sum_a_data_need_filter1=sum_a_data_need;
% sum_a_data_need_filter1(:,1:500)=0;
sum_a_data_need_filter=sum_a_data_need_filter1;
% [m,n]=size(sum_a_data_need_filter1);
%%
sum_a_data_need_filter_move=zeros(a,length(sum_a_data_need(1,:)));
for i=1:a
    for k=1:length(sum_a_data_need(1,:))
    sum_a_data_need_filter_move(i,round(k/2))=sum_a_data_need_filter(i,k);
    end
end
% %% 信号平移
% sum_a_data_need_filter_move=circshift(sum_a_data_need_filter,[1 -50]);
%%
sum_a_data_need_filter_move(sum_a_data_need_filter_move==-Inf)=0;
sum_a_data_need_filter_move(sum_a_data_need_filter_move==0)=min(min(sum_a_data_need_filter_move));
sum_a_data_need_filter_move=sum_a_data_need_filter_move-min(min(sum_a_data_need_filter_move));
%% 通过窗的大小调节滤波的范围
sum_new_data=zeros(256,3200);
win=1;%选4
for j=1:256
    new_data=zeros(1,3200);
    aa=sum_a_data_need_filter_move(j,:);
    for i=1:length(aa)/win
        loc=aa(1+win*(i-1):win*i);
        loc_max=find(loc==max(loc));
        new_data(loc_max+(i-1)*win)=max(loc);
    end
    sum_new_data(j,:)= new_data;
end 
sum_new_data(:,1400:1700)=0;
% %% 
% sum_a_data_need_end=zeros(a,length(sum_a_data_need(1,:)));
% for kk=1:a
%     aa=sum_a_data_need_filter_move_hilbert(kk,:);
%     [pks,locs] = findpeaks(aa); 
%     sum_a_data_need_end(kk,locs)=pks;
% end
% sum_a_data_need_end=circshift(sum_a_data_need_end,[1 -40]);%%20-30
sum_new_data=circshift(sum_new_data,[1 -40]);%%20-30
sum_new_data(:,1400:3200)=0;
%% 包络
[m,n]=size(sum_new_data);
for i=1:m
    x=sum_new_data(i,:);
    t_end=307e-6;
    t=t_end/n:t_end/n:t_end; 
    z=hilbert(x');                          % 希尔伯特变换
    baoluo=abs(z);                               % 包络线
    sum_new_data_hilbert(i,:)= baoluo';
end
figure(1)
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1), pos(2)-100,pos(3),pos(4)]);
plot(t,x,'k'); hold on;
plot(t,a,'r--','linewidth',2);
title('包络线'); ylabel('幅值'); xlabel(['时间/s' 10 '(a)']);
% ylim([-2,2]);
% mesh(sum_new_data);
sum_new_data_hilbert(:,1400:3200)=0;
image(sum_new_data_hilbert);
axis square
%% 移动平均滤波器
% t = 1:3200;
% % rng default  %initialize random number generator
% % x = sin(t) + 0.25*rand(size(t));
% windowSize = 20; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% for i=1:256
%     y = filter(b,a,sum_new_data_hilbert(i,:));
% end
% %%
% plot(t,sum_new_data_hilbert(1,:))
% hold on
% plot(t,y)
% legend('Input Data','Filtered Data')
