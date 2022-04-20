clear,clc%,close all

filename = "adc_1280Hz_1.txt";
%filename = "adc.txt";
Fd=1280;% Частота дискретизации (Гц)
%Fd / DIV = 300 Hz -> 150 Hz - частота среза (Найквиста)
FftL=2048;% Количество линий Фурье спектра
N_all=1536;
T = 1; %sec
N = T*Fd; %80;
DIV = 16;

fileID = fopen(filename,'r');
adc1 = fscanf(fileID, '%d', N_all)';
adc2 = fscanf(fileID, '%d', N_all)';
adc3 = fscanf(fileID, '%d', N_all)';
adc4 = fscanf(fileID, '%d', N_all)';
fclose(fileID);


adc1 = adc1 - adc2;
adc4 = adc3 - adc4;
adc2 = zeros(1,N_all);
adc3 = zeros(1,N_all); 


t = 0:1/Fd:(N_all-1)/Fd;

%% test
if 0
    Fs = 10;
    A1 = 2;
    adc1 = 10*cos(2*pi*Fs*t);% + A1*cos(pi*Fs*t) - A1*sin(pi*Fs*t); 
    adc4 = 10*sin(2*pi*Fs*t);% - A1*sin(pi*Fs*t); 
end


figure(20)
plot(t,adc1,t,adc2,t,adc3,t,adc4)
title("до фильтрации")

%% FIR-filter
%[adc1, adc4] = fir_filter(adc1,adc4,N, Fd);
fir_length = 9;
h_low = fir1(fir_length-1, 0.05,hamming(fir_length));
%h = [-0.00197114084026180210948431259510016389 -0.005329221856815129568230027246045210632 -0.013471051671775362279515064756196807139 -0.021652088787228292859898814981534087565 0.976921763735987380705694249627413228154 -0.021652088787228292859898814981534087565 -0.013471051671775362279515064756196807139 -0.005329221856815129568230027246045210632 -0.00197114084026180210948431259510016389]; 

h = [   %частота разворота 80 Гц, среза 1 Гц
-0.00340725044521170746664173734075120592 
-0.006951756519563004506345738064965189551
-0.018435911117199745706818347912303579506
-0.042278334890664091838541338574941619299
-0.087523822861556813124117581992322811857
-0.185762852473975587086485461441043298692
-0.625625959460234359887920163600938394666
0.625625959460234359887920163600938394666
0.185762852473975587086485461441043298692
0.087523822861556813124117581992322811857
0.042278334890664091838541338574941619299
0.018435911117199745706818347912303579506
0.006951756519563004506345738064965189551
0.00340725044521170746664173734075120592 
];
if Fd == 80
    adc1 = conv(adc1, h, 'valid');
    adc4 = conv(adc4, h, 'valid');
end

%% среднее значение, чтобы убрать постоянную составляющую
%adc_M = [-101.708984375000 -43.1269531250000 -29.5234375000000 -111.992187500000]; %среднее на 40Гц
adc_M = [-100.1640625000000 -43.066406250000000 -27.430664062500000 -111.1464843750000];% среднее без магнитов
%adc_M = [39.4414062500000 0 0 -29.4036458333333];   %diff mode
filter_type = 3;
if Fd == 1280
    filter_type = 2;
end

if filter_type == 0    %вычесть измеренное значение без магнита, либо среднее из N_all значений
    adc1 = adc1 - adc_M(1);
    adc2 = adc2 - adc_M(2);
    adc3 = adc3 - adc_M(3);
    adc4 = adc4 - adc_M(4);  
elseif filter_type == 1
    [adc1,adc_M(1)] = normalize_adc(adc1);
    [adc2,adc_M(2)] = normalize_adc(adc2);
    [adc3,adc_M(3)] = normalize_adc(adc3);
    [adc4,adc_M(4)] = normalize_adc(adc4);  
end

% выкинуть часть значений
adc1 = adc1(1:N);
adc2 = adc2(1:N);
adc3 = adc3(1:N);
adc4 = adc4(1:N);
t = 0:1/Fd:(N-1)/Fd;


%% прореживание + накопление
adc1_f = zeros(1,N/DIV);
adc2_f = zeros(1,N/DIV);
adc3_f = zeros(1,N/DIV);
adc4_f = zeros(1,N/DIV);
adc1_1 = zeros(1,N/DIV);
adc2_1 = zeros(1,N/DIV);
adc3_1 = zeros(1,N/DIV);
adc4_1 = zeros(1,N/DIV);
j = 1;
for i=1:DIV:N
    adc1_1(j) = adc1(i);
    adc2_1(j) = adc2(i);
    adc3_1(j) = adc3(i);
    adc4_1(j) = adc4(i);
    for k = 0:DIV-1
        adc1_f(j) = adc1_f(j) + adc1(i+k);
        adc2_f(j) = adc2_f(j) + adc2(i+k);
        adc3_f(j) = adc3_f(j) + adc3(i+k);
        adc4_f(j) = adc4_f(j) + adc4(i+k);
    end
    j = j + 1;
end

if filter_type == 2
    adc1 = conv(adc1, h, 'valid');
    adc4 = conv(adc4, h, 'valid');
    adc1_1 = conv(adc1_1, h, 'valid');
    adc4_1 = conv(adc4_1, h, 'valid');
    adc1_f = conv(adc1_f, h, 'valid');
    adc4_f = conv(adc4_f, h, 'valid');
end

t = 0:1/Fd:(length(adc1)-1)/Fd;
%[adc1_f,adc2_f] = fir_filter(adc1_f,adc2_f,N/DIV, Fd/DIV);
%[adc3_f,adc4_f] = fir_filter(adc3_f,adc4_f,N/DIV, Fd/DIV);

%% comlex signal
adc_h0 = adc4 + 1i*adc1;
adc_h1 = adc4_1 + 1i*adc1_1;
adc_h2 = adc4_f + 1i*adc1_f; 


phi0_wrap = angle(adc_h0);
phi1_wrap = angle(adc_h1);
phi2_wrap = angle(adc_h2);
phi0 = unwrap(angle(adc_h0));
phi1 = unwrap(angle(adc_h1));
phi2 = unwrap(angle(adc_h2));

%% Мгновенная частота
freq0 = freq_average(phi0, Fd);
freq1 = freq_average(phi1, Fd/DIV);
freq2 = freq_average(phi2, Fd/DIV);



%% Оценка Сигнал/шум
%S_noise = S_signal без вращения крыльчатки
%S_noise = [3.232887875554642e+04 8.326142276274370e+03 3.739248616969993e+04];
S_noise = [1.295231892604766e+05 3.149973302814137e+04 1.275036363603583e+05];  %диф. режим 1280Гц + ФВЧ    27,15 dB с накоплением
%S_noise = [2.250310866490334e+04 5.667883942196749e+03 2.933996618414966e+04];  %диф. режим 1280Гц ФВЧ и ФНЧ   32,16 dB с накоплением
%S_noise = [9.073885273070184e+04 2.262090713164671e+04 7.979012391432036e+04];  %single 1280 + ФВЧ   25,73 dB с накоплением
S_signal = [0 0 0];
snr_dB = [0 0 0];
%S_noise(1) = 3.130831005069404e+04;
%S_noise(2) = 8.490421440784590e+03;   % прореживание
%S_noise(3) = 2.811628502007019e+04;   % накопление
[snr_dB(1),S_signal(1)] = snr_value(adc_h0,S_noise(1),Fd,FftL);
[snr_dB(2),S_signal(2)] = snr_value(adc_h1,S_noise(2),Fd/DIV,FftL);
[snr_dB(3),S_signal(3)] = snr_value(adc_h2,S_noise(3),Fd/DIV,FftL);
psd1_value(adc_h2,Fd/DIV);
%psd1_value(adc1_f,Fd/DIV);

 
%% Graphics
set(0,'defaultAxesXLimSpec', 'tight')
set(0,'defaultAxesYLimSpec', 'tight')

figure(11)
subplot(3,1,1);
plot(phi0_wrap, '-o')
title('Фаза');
subplot(3,1,2);
plot(phi1_wrap, '-o')
subplot(3,1,3);
plot(phi2_wrap, '-o')
figure(12)
subplot(3,1,1);
plot(phi0, '-o')
title('Фаза unwrap');
subplot(3,1,2);
plot(phi1, '-o')
subplot(3,1,3);
plot(phi2, '-o')

figure(1)
subplot(3,1,1);
plot(t,adc1,'-o',t,adc4,'-o')
title("Исходный сигнал " + freq0 + " Гц, SNR " + snr_dB(1) + " dB");
subplot(3,1,2);
t1 = t(1:length(adc1_1));
plot(t1,adc1_1,'-o',t1,adc4_1,'-o')
title("прореживание " + freq1 + " Гц, SNR " + snr_dB(2) + " dB"); 
%{
 plot(F,A1,F,A2)
title('преобразование Фурье'); 
%}
subplot(3,1,3);
t2 = t(1:length(adc1_f));
plot(t2,adc1_f,'-o',t2,adc4_f,'-o')
title("накопление " + freq2 + " Гц, SNR " + snr_dB(3) + " dB"); 


%% FIR-filter с прореживанием в 2 раза
function [i_f,q_f] = fir_filter(inph,quadr,N, Fd) 

    N_filter = 11;
    h_filter = zeros(1,N_filter);
    i_f = zeros(1,(N-N_filter+1)/2);
    q_f = zeros(1,(N-N_filter+1)/2);
    centre = (N_filter+1)/2;

    f_cut = 200;
    %h = fir1(N_filter-1, f_cut/(Fd/2));

    h(1) = -0.001951997961388574889460278960484629351;%-0.000000000000000001551078847964774085357;
    h(2) =  0.013377862194433463166598485827307740692;%-0.022663985459552640072677931470934709068;
    h(3) = -0.016179305194684238944358156686575966887;%0.000000000000000010469782223762223728009;
    h(4) = -0.076887053313202766147149702646856894717;%0.273977082565523999413414912851294502616;
    h(5) =  0.268658314486578342350497905499651096761;%0.497373805788057288257419941146508790553;
    h(6) =  0.625964359576527740181006720376899465919;%0.273977082565523999413414912851294502616;
    h(7) =  0.268658314486578342350497905499651096761;%0.000000000000000010469782223762223728009;
    h(8) = -0.076887053313202766147149702646856894717;%-0.022663985459552640072677931470934709068;
    h(9) = -0.016179305194684238944358156686575966887;%-0.000000000000000001551078847964774085357; 
    h(10) =  0.013377862194433463166598485827307740692;
    h(11) = -0.001951997961388574889460278960484629351;

    i_f = conv(inph, h, 'valid');
    q_f = conv(quadr, h, 'valid');

    %{
    for i = 0:2:(N-N_filter)
        i_f(i/2+1) = inph(i+centre)*h(centre);
        for j = 2:2:N_filter
            q_f(i/2+1) = q_f(i/2+1) + quadr(i+j)*h(N_filter+1-j);
        end
    end 
    %}

end


%% Мгновенная частота
% количество точек на промежутке 2*pi
% delta_phi - N-1 точек
% 2*pi      - N1
% f = phi / (2*pi*dt*N) = 1 / (dt * N1)
function freq1 = freq_average(phi, Fd)

    %{
        sum1 = sum(phi0);
    sum1 = sum1 * 2 / (length(phi0)-1);
    sum1 = sum1 / 2 / pi;

    fr = sum1 * Fd / N;
    %}
    N = length(phi);  %количество интервалов
    if 1
        sum_phi = sum(phi) / N / pi;
        freq1 = sum_phi * Fd / (N - 1);
    else
        N1 = 2*pi*(N-1)/(phi(end) - phi(1));
        freq1 = Fd / N1;
    end
end


%% Мгновенная частота - производная мгновенной фазы
% w = dphi/dt = (phi(n+1)-phi(n))*Fd = 2*pi*f	(w - угловая частота, f - мгновенная частота)
function freq = freq_calc(phi, Fd)
    
    freq = zeros(1, length(phi));
    for m=1:length(phi)-1
        freq(m) = Fd/(2*pi)*(phi(m+1)-phi(m));
    end
    freq(end) = freq(end-1);
end

%% Оценка Сигнал/шум
function [w,psd_dB] = psd1_value(adc, Fd)
    %snr(adc, Fd, 6);
    %periodogram(adc);  %оценка спектральной плотности мощности, основанная на вычислении квадрата модуля преобразования Фурье
    [p_adc, w] = periodogram(adc,[],[],Fd); % compute the psd(power spectral density) estimate
    np_adc = p_adc / max(p_adc); % normalize
    psd_dB = 10 * log10(np_adc);
    
    figure(100)
    plot(w, psd_dB)
    grid on
    xlabel('Hz')
    ylabel('dB')
end


%% Оценка Сигнал/шум
function [snr_dB,S_signal] = snr_value(adc,S_noise,Fd,FftL)
    adc = [adc, zeros(1,FftL-length(adc))];
    A = abs(fftshift(fft(adc, FftL)));
    f=-Fd/2:Fd/FftL:Fd/2-1/FftL;% Массив частот вычисляемого спектра Фурье
    %f = [-FftL/2+1:FftL/2];
    figure(23)
    plot(f,A);
    xlabel('Hz')

    S_signal = sum(A); %площадь АЧХ сигнала 

    snr_dB = 20*log10(S_signal/S_noise);
end

%% FFT
function [F,A] = fft_norm(adc,Fd,FftL)
    
    F=0:Fd/FftL:Fd/2-1/FftL;% Массив частот вычисляемого спектра Фурье
    A = abs(fft(adc, FftL));
    A = A(1:length(F));
    %A1=2*A1./FftL;% Нормировка спектра по амплитуде
    %A1(1)=A1(1)/2;% Нормировка постоянной составляющей в спектре
    figure(101)
    plot(F,A)
    title('преобразование Фурье'); 
end

%% вычитаем среднее значение, чтобы убрать постоянную составляющую
function [adc,adc_M] = normalize_adc(adc)
    adc_M = mean(adc);
    adc = adc - adc_M;
end
