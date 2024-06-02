clear
close all
clc
c=1540;

fs=25e6;     
N_elements=128;
img_width=38/1000;
d_x=img_width/N_elements; %  Increment for image（每一个振元有多长）
min_sample=0;

load RF_Filter.mat
data1=RF_Filter;%导入射频数据

for y=1:128
	rf_env=abs(hilbert(data1(:,y)));
	rf_env=rf_env';
	env(1:max(size(rf_env)),y)=rf_env;
end

D=3;
dB_range=40;  % Dynamic range for display in dB
disp('Finding the envelope')

log_env=env(1:D:max(size(env)),:)/max(max(env));  
log_env=log(log_env+0.01);%%%0.6复合：0.0003，单帧：0.003。 0:复合0.001，单帧：0.01 compensation:0.001(复合)
log_env=log_env-min(min(log_env));
log_env=80*log_env/max(max(log_env));

[n,m]=size(log_env);
ID_bmode=round(n/N_elements); 
new_env=zeros(n,m*ID_bmode);
for u=1:n
    if(rem(u,100) == 0)
    end
    new_env(u,:)=abs(interp(log_env(u,:),ID_bmode));
end

[n,m]=size(new_env);
fn_bmode=fs/D;

figure(1)
h=image(((1:(ID_bmode*N_elements-1))*d_x/ID_bmode-N_elements*d_x/2)*1000,((1:n)/fn_bmode+min_sample/fs)*1540/2*1000,new_env);
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(gray)
ylim([0 38])
%% 计算FWHM
[n,m]=size(new_env);
point_coordinate = [0 0 60];
phntPoints=point_coordinate/1000;
img_width=38/1000;
maxSamples=n;
nBeams=m;

img=new_env;

img = img / max(img(:));
img = log(img+0.01);
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
figure(1)
plot(scaled_data,img(maxRow,:));%%(418:602)
xlabel('x [mm]'); ylabel('PSF/PSFmax');%
%  axis([(phntPoints(1,1)-imageXDim/2)*10000 (phntPoints(1,1)+imageXDim/2)*10000 0 1.05]);
title(['PSF(Lateral Profile); 横向FWHM = ' num2str(FWHM_lateral)]);
% grid on
