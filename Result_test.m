%%
clc
clear
close all
load sensor_data.mat
%% 幅度谱
figure
data = zeros(1664,128);
for k=1:128
    data(:,k)=sensor_data(:,k,1);
end
% subplot(2,1,1)
SpecPlot(data,25e6,'g',2);
hold on
%%
SpecPlot(RF(:,1),25e6,'b',2);
hold on
%%
SpecPlot(RF(:,1),25e6,'r',2);
legend('汉宁窗合成','DOA合成后')
%% 时域谱
figure
% subplot(2,1,2)
x=hilbert(data);
plot(1:H,abs(x),'g');
hold on
legend('合成前')
y=hilbert(RF(:,1));
plot(1:H,abs(y),'r');
xlabel('采样点');
ylabel('幅度');
%% 图像分析
clc
% A = imread('img.png'); %换成自己的图
%load new_env.mat
A=new_env;
GRAY_A = double(A).*255;
ENG_GRAY_A=GRAY_A.*GRAY_A;%计算能量
[x1,y1] = size(ENG_GRAY_A);
X = 0:x1-1;
Y = 0:y1-1;

if(size(ENG_GRAY_A,3)>1)

ENG_GRAY_A=rgb2gray(ENG_GRAY_A);

end

ENG_GRAY_A=double(ENG_GRAY_A);
figure
subplot(2,1,1)
mesh(flipdim(ENG_GRAY_A,1)),title('灰度能量图')
subplot(2,1,2)
mesh(flipdim(ENG_GRAY_A,1)),title('灰度能量图')

%% 保存图窗图像
f=getframe(gca);
imwrite(f.cdata,'img.png')
%%
figure
data = zeros(1664,128);
for k=1:128
    data(:,k)=sensor_data(:,k,1);
end

x=hilbert(RF(:,1));
plot(1:H,abs(x),'g');
hold on
legend('合成前')
%%
y=hilbert(RF(:,1));
plot(1:H,abs(y),'r');
xlabel('采样点');
ylabel('幅度');
hold on
rectangle('Position',[0 0 200 3000],'EdgeColor','b')