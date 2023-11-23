% clear all;clc
% %select a file,import it and display it as an image
% [filename filepath]=uigetfile({'*.bmp','bmp file(*.bmp)';'*.*','All files(*.*)'},'Please select a file');
% % [filename filepath]=uigetfile({'*.jpg','bmp file(*.bmp)';'*.*','All files(*.*)'},'Please select a file');
% if isequal(filename,0)
%    msgbox('User selected Cancel','Warning','warn');
%    return
% else if ~isequal(filename(1,length(filename)-3:length(filename)),'.bmp')
%         msgbox('Wrong file type','Warning','warn');
%         return
%     else
%        msgbox(['User selected ',[filepath, filename]])
%     end
% end
% A=imread([filepath filename],'bmp');
% clear filename filepath
% function A = myfun(scan_lines_fund,[])
% load G:\1\scan_lines_fund1
% RFData_256_SA1(1:70,:)=0;
A=scan_lines_fund;
A(:,1:20)=0;
A(:,1)=max(A(:));
% A(1:500,290:340)=A(1:500,230:280);
% A(501:1000,330:390)=A(501:1000,230:290);
% A(1001:1500,340:400)=A(1001:1500,230:290);
% A(1501:1700,330:380)=A(1501:1700,230:280);
% A(1701:1800,320:360)=A(1701:1800,230:270);
% A(1801:2047,300:340)=A(1801:2047,230:270);
% A(128:130,:)=0;
% A(1950:1970,:)=A(1930:1950,:);
% A(2020:2047,:)=A(21:48,:);%%%timo
B=A(:,1:1600); 
A=rot90(B,2);
figure;
imshow(A,[]);

rM=1000; %the rectangle image is of the size rM*rM
B=zeros(rM);

[M,N]=size(A);
deltR=N/((rM-1)/2);%矩形图上没点所代表的半径长度相当于原图上deltR个点代表的长度
delta=2*pi/M;

Ox=(rM-1)/2;%the position of the orgin of the rectangle image
Oy=(rM-1)/2;
for i=1:1:rM
    for j=1:1:rM
        rx=j-Ox;ry=i-Oy;
        r=sqrt(rx^2+ry^2);
        if r<=(rM-1)/2
            rmin=floor((r-1)*deltR+1);
            rmax=ceil(r*deltR);
            if rmin<1  %matrix indices couldnot be out of the range
                rmin=1;
            end
            if rmax>N
                rmax=N;
            end
        
            angle=atan2(ry,rx); 
            if angle<0
                angle=angle+2*pi;
            end
            amin=floor(angle/delta);
            amax=ceil(angle/delta);
            if amin<1
                amin=1;
            end
            if amax>M
                amax=M
            end
            
            %找到新图上每点在原图上对应的区域后，对该区域的值求和然后求平均以得
            %到新图上对应点的值
            sum=0;count=0;
            for m=rmin:1:rmax
                for n=amin:1:amax
                    sum=double(A(n,m))+sum;
                    count=count+1;
                end
            end
            if count~=0
                B(i,j)=round(sum/count);
    
            end
        end
    end
end
% clear M N Ox Oy i j rx ry r rmin rmax amin amax angle sum count
% clear deltR delta m n
% B(499:503,500:100)=B(492:497,500:100);
%%
% env_dB=30*log10(B);
% env_dB=env_dB-max(max(env_dB));
% B=127*(env_dB+60)/200;
%%
figure
imshow(B,[])
figure(3)
ll=length(B(1,:));
x=[200/ll:200/ll:200]-100;
imagesc(x,x,B,[0,50]);%20
axis square
colormap(gray)
xlabel('y-position[mm]');
ylabel('x-position[mm]');
% hleg1 = legend('Sampling rate 20%','Sampling rate 40%','Sampling rate 60%','Sampling rate 80%','Orientation','vertical','box','off');%vertical竖直的、horizontal水平的
% set(hleg1, 'Position', [0.629367562364342 0.683921396461717 0.217187494318932 0.176195421114781]);
set(gca,'FontName','Times New Roman','fontweight','bold','FontSize',16,'LineWidth',1);