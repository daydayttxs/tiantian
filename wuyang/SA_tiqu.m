clear all;clc;%%??3?
load D:\实验室\Bmodule\实测脚本\4.11\RFdata_11_Apri_shengsu_5.3_SA20.mat
RcvData1=RcvData{1}(:,:,1);
% load E:\user-MATLAB\Test\实测\256合成孔径间隔1\data\原始数据\RcvData_water_4
% RcvData1=RcvData_water_4{1};
% data_need_1=data_need(1:3200,:);
% data_need_1=double(data_need_1);
RFData_yuan=double(RcvData1);  
nSample=3200;
nEl=256;
SA=1;
RFData=zeros(nSample,nEl);
for i=1:nEl
    data=zeros(nSample,nEl);
    if i+SA-1<=256
        data(:,i:i+SA-1)= RcvData1(1+(i-1)*nSample:nSample+(i-1)*nSample,i:i+SA-1);
    else
        data(:,i:nEl)=RcvData1(1+(i-1)*nSample:nSample+(i-1)*nSample,i:nEl);
        data(:,1:i+SA-1-nEl)=RcvData1(1+(i-1)*nSample:nSample+(i-1)*nSample,1:i+SA-1-nEl);
    end
    RFData=RFData+data;
end  
RFData_SA1=RFData;
% RFData_sig(1:500,:)=0;
mesh(20*log(abs(RFData)));
view(0,90)