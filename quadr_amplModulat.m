clear,clc,close all

%% Параметры
N=64;% Количество точек дискретизации
Fs=100;% Частота синусоиды (Гц)
Fd=512;% Частота дискретизации (Гц)
FftL=512;
Tm=N/Fd;% Длина сигнала (с)
t=0:1/Fd:(Tm-1/Fd);% Массив отсчетов времени
fi_1 = 0;%pi/6;

inph=cos(2*pi*Fs*t+fi_1) + 0.2*cos(pi*Fs*t) - 0.2*sin(pi*Fs*t);
quadr=sin(2*pi*Fs*t+fi_1) - 0.2*sin(pi*Fs*t);

%% Noise = 0
A_noise = 0.1;
inph = inph + A_noise*randn(1,N);
quadr = quadr + A_noise*randn(1,N);

%% Phase
phi1 = angle(inph + 1i*quadr);
phi = unwrap(phi1);

%% Мгновенная частота
N1 = 2*pi*(N-1)/(phi(end) - phi(1));
freq1 = Fd / N1;


%% Graphics
% FFT
F=0:Fd/FftL:Fd/2-1/FftL;% Массив частот вычисляемого спектра Фурье
s_fft = abs(fft(inph, FftL));
s_fft = s_fft(1:length(F));
plot(F,s_fft)
title('преобразование Фурье');

figure()
subplot(3,1,1)
plot(t,inph,t,quadr)
title("Мгновенная частота " + freq1 + " Гц (реальная " + Fs + "Гц)");
subplot(3,1,2)
plot(t,phi1)
title('Фаза');
subplot(3,1,3)
plot(t,phi)
