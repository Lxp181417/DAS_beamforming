function SpecPlot(signal,Fs,color,choose)
% SpecPlot 绘制信号单边频谱
% signal：输入信号
% Fs：信号的采样频率
% color：曲线颜色
% choose：绘制模式选择
% choose=1：绘制原始曲线；choose=2：绘制平滑曲线；choose=3：同时绘制两者；

% FFT 变换
N = length(signal);              % 信号长度
X = fft(signal);                 % FFT 变换
X_mag = abs(X);                  % 幅度谱

% 计算频率轴
f = (0:N-1)*(Fs/N);              % 频率轴

% 计算单边幅度谱
X_mag_single_side = X_mag(1:N/2+1);  % 保留单边幅度谱部分

% 绘制幅度谱
n=2;
if choose==1
    plot(f(1:N/2+1), X_mag_single_side,color);
elseif choose==2
    for i=1:n
        X_mag_single_side=smooth(X_mag_single_side);%数据平滑
    end
    plot(f(1:N/2+1), X_mag_single_side,color);
elseif choose==3
    plot(f(1:N/2+1), X_mag_single_side,'b');
    hold on

    for i=1:n
        X_mag_single_side=smooth(X_mag_single_side);%数据平滑
    end
    plot(f(1:N/2+1), X_mag_single_side,color);
end

xlabel('频率 (Hz)');
ylabel('幅度');
end