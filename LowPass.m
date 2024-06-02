function Hd = LowPass(Fc,Fs)
%LOWPASS 返回离散时间滤波器对象。

% MATLAB Code
% Generated by MATLAB(R) 9.11 and Signal Processing Toolbox 8.7.
% Generated on: 23-Apr-2024 17:29:53

% FIR Window Lowpass filter designed using the FIR1 function.

% All frequency values are in kHz.
% Fs: Sampling Frequency
% Fc: Cutoff Frequency

N    = 10;       % Order
flag = 'scale';  % Sampling Flag

% Create the window vector for the design algorithm.
win = rectwin(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);

% [EOF]