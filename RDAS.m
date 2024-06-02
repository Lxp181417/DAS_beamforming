%% 获取通道数据
clc
clear
close all
data=load("JDMRcvData.mat");
% fieldnames(data)%查看结构体元素
RcvData=getfield(data, 'RcvData');%获取结构体中的数据

%数据整理
sensor_data = zeros(1664,128,128);
for e=1:128
    for k=1:128
        sensor_data(:,e,k)=RcvData(1+1664*(k-1):1664+1664*(k-1),e);
    end
end
save sensor_data.mat
%% 导入通道数据
clc
clear
close all
load sensor_data.mat
%% 参数定义
% 采样频率：Fs=25e6 
% 阵元间距：pitch=30*1e-6 
% 声速：    c=1540m/s 
% 发射阵元：N=128
clc
% close all
clearvars -except sensor_data %清除变量
N=128;            % 超声换能器阵列单元数量
W=N;              % 定义图像大小
H=1664;
RF=zeros(W,H);    % 定义射频信号大小
DOA_window=zeros(H,W);
RF_Filter=zeros(H,W);

fs=25e6;          % 数据的采集频率，单位：Hz
f0=6.25e6;        % 换能器中心频率，单位：Hz
dt=1/fs;          % 数据的采集时间分辨率，单位：s
c=1540;           % 介质声速，单位：m/s

pitch=30e-6;      % 阵列单元间隔，单位：m
w=pitch;          % 定义图像像素尺寸，单位：m
h=pitch;      % 定义图像像素尺寸，单位：m

lambda=c/6.25e6;  % 定义波长

%% 计算窗函数
for j=1:H
    for k=1:N
        len=N/2*w-k*w;
        depth=j*h;
        theta=atan(len/depth);
        DOA_window(j,k)=(cos(2*pi*pitch*sin(theta)/lambda)+1)*cos(theta)/2;
    end
end
Lateral_window=hanning(W);%阵元方向加窗--横向
Axial_window=hamming(H);%探测深度方向加窗--轴向

%% RDAS波束合成
tic%计时
for i=1:W
    for j=1:N
        for k=1:H
            index=(c/(fs*2*h))*k-(fs*w^2/(c*2*h))*(i-j)^2/k;
            index=round(index);
            if index<=0
                index=1;
            elseif index>H
                index=H;
            end

            if(index<=size(sensor_data,1))
                %DOA
                RF(i,index) = RF(i,index) + sensor_data(k,j,i).*Lateral_window(j).*Axial_window(k).*DOA_window(k,j);
%                 RF(i,index) = RF(i,index) + sensor_data(k,j,i).*Lateral_window(j);
            end
        end
    end
end
RF=RF'; % 图像旋转

toc
disp(['运行时间: ',num2str(toc)]);

%% 加窗/滤波
cutoff_freq=9.8e6;  % 高通滤波器截止频率，单位：Hz

for i=1:W
    RF_Filter(:,i)=filter(HighPass(cutoff_freq,fs),RF(:,i));% 高通滤波
%     RF_Filter(:,i)=filter(LowPass(cutoff_freq,fs),RF(:,i));% 低通滤波
%     RF_Filter(:,i)=filter(BandPass(),RF(:,i));% 带通滤波
end

%% 成像 
% Compress the data to show 60 dB of dynamic range for the cyst phantom
% image
%  Clibrated 60 dB display made
f0=6.25e6; % Transducer center frequency [Hz]，超声换能器的中心频率
fs=25e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
no_lines=128; % Number of lines in image，图像中的扫描线数量
image_width=38/1000; % Size of image sector，图像扇区的大小
% d_x=image_width/no_lines; % Increment for image，图像增量
N_elements=128; % 阵元的数量
d_x=image_width/N_elements; % Increment for image，图像的增量
% Read the data and adjust it in time，读取数据并调整时间
min_sample=0;

rf_data = RF_Filter; % 导入射频信号
% rf_data = RF;

for i=1:no_lines
    % Find the envelope，找到包络
    rf_env=abs(hilbert(rf_data(:,i))); % 使用希尔伯特变换获取包络
    env(1:max(size(rf_env)),i)=rf_env;
end

% Do logarithmic compression，进行对数压缩
D=3; % Compression factor，压缩因子
log_env=env(1:D:max(size(env)),:)/max(max(env)); % 对包络进行压缩
log_env=log(log_env+0.3); % 取对数并添加常数以避免出现零值
log_env=log_env-min(min(log_env)); % 最小值规范化
log_env=80*log_env/max(max(log_env)); % 设定动态范围为60dB，将对数值映射到0-80之间

G= im2double (log_env); % Convert the image to double precision，将图像转换为双精度

% Make an interpolated image，生成插值图像
[n,m]=size(log_env);
ID_bmode=round(n/N_elements); % Compute interpolation factor，计算插值因子
new_env=zeros(n,m*ID_bmode); % 创建新的图像矩阵

for u=1:n
    if(rem(u,100) == 0) % 如果u能被100整除
    end
    new_env(u,:)=abs(interp(log_env(u,:),ID_bmode)); % 对每行进行插值
end

[n,m]=size(new_env);
% new_env1=new_env(250:505,:);
fn_bmode=fs/D; % Compute sampling frequency for display，计算显示的采样频率
% figure(1) clf
h=image(((1:(ID_bmode*N_elements-1))*d_x/ID_bmode-N_elements*d_x/2)*1000,((1:n)/fn_bmode+min_sample/fs)*1540/2*1000,new_env); % Display the image，显示图像
colormap(gray); % Set the colormap to grayscale，设置色图为灰度
brighten(+0.3) % Increase the brightness of the image，增加图像的亮度
ylim([0 38])
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
% saveas( 1, 'image.png'); %保存Figure1窗口的图像

%% 计算FWHM
[n,m]=size(new_env);
point_coordinate = [0 0 60];
phntPoints=point_coordinate/1000;
img_width=38/1000;
maxSamples=n;
nBeams=m;

img=new_env;

img = img / max(img(:));
img = log(img + 0.01);
img = img - min(img(:));
img = img / max(img(:));
imgXDim=img_width/nBeams;
dxBeamCoord=d_x;
%-- 查找图像最大值<maxPixel>所在的列<maxColumn>
maxPixel=0;
maxRow=0;
for i=1:maxSamples
    for j=1:nBeams
        if maxPixel<img(i,j)
            maxPixel=img(i,j);
            maxRow=i;
            maxColumn=j;
        end
    end
end
%-- 查找列<maxColumn>轴向方向上第一个非零值<nonzeroRow1>
for i=maxRow:-1:1
    if img(i,maxColumn)>0
        nonzeroRow1=i;
    end
end
%-- 查找列<maxColumn>轴向方向上（从信号的另一侧）的非零值<nonzeroRow2>
for i=maxRow:maxSamples
    if img(i,maxColumn)>0
        nonzeroRow2=i;
    end
end
%-- 轴向/径向FWHM
axial_res = zeros(maxSamples, 1);
for i=2:(nonzeroRow1-1)
    axial_res(i,1)=axial_res((i-1),1)+img(nonzeroRow1,maxColumn)/(nonzeroRow1-1);
end
axial_res(nonzeroRow1:nonzeroRow2,1)=img(nonzeroRow1:nonzeroRow2,maxColumn);
for i=(nonzeroRow2+1):maxSamples
    axial_res(i,1)=axial_res((i-1),1)-img(nonzeroRow2,maxColumn)/(maxSamples-(nonzeroRow2+1));
end
cas_rezu=0;
pridavek_na_ose_z=0;
while cas_rezu<0.5
    cas_rezu = cas_rezu + ((1/fs) * c)/2* 1000;
    pridavek_na_ose_z = pridavek_na_ose_z+1;
end
%-- 在轴向上查找大于图像最大值<maxPixel>一半的第一个值<axial_sub1>
i=maxRow;
axialTemp1 = maxPixel;
while axialTemp1 >= (maxPixel/2)
    axialTemp1 = img(i,maxColumn);
    axialSub1 = i;
    i=i-1;
end
%-- 查找第二个值<axial_sub2>
i=maxRow;
axialTemp2 = maxPixel;
while axialTemp2 >= (maxPixel/2)
    axialTemp2 = img(i,maxColumn);
    axialSub2 = i;
    i=i+1;
end

%-- 搜索列号
for i=1:(nBeams/2)
    if img(maxRow,i)<(maxPixel/2)
        lateralSub1=i;
    end
end
for i=nBeams:-1:(nBeams/2)
    if img(maxRow,i)<(maxPixel/2)
        lateralSub2=i;
    end
end
xAxis=((phntPoints(1,1)-(imgXDim/2))+[0:nBeams-1]*dxBeamCoord) * 1000;
zAxis = (([0:maxSamples-1] / fs) * c)/2* 1000;
[m1,n1]=size(xAxis);
[m2,n2]=size(zAxis);
FWHM_lateral=abs(lateralSub1-lateralSub2)*38/n1;
FWHM_axial=abs(axialSub2-axialSub1)*max(zAxis)/n2;

side_lobe_height1=mean(img(nonzeroRow1:nonzeroRow2,1:lateralSub1-1),'all');
side_lobe_height2=mean(img(nonzeroRow1:nonzeroRow2,lateralSub1+1:end),'all');
side_lobe_height=(side_lobe_height1+side_lobe_height2)/2;
disp(['横向FWHM:' num2str(FWHM_lateral)]);
disp(['旁瓣高度:' num2str(side_lobe_height)]);


% 原始数据范围
min_original = min(xAxis);
max_original = max(xAxis);

% 目标数据范围
min_target = -19;
max_target = 19;

% 线性变换
scaled_data = ((xAxis - min_original) / (max_original - min_original)) ...
              * (max_target - min_target) + min_target;

% Plot
plot(scaled_data,img(maxRow,:));%%(418:602)
xlabel('x [mm]'); ylabel('PSF/PSFmax');%
title(['PSF(Lateral Profile); 横向FWHM = ' num2str(FWHM_lateral)]);
% grid on
