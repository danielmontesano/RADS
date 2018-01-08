function [ y_dem ] = demoduladorExpansion( y, Rb,f0,f1)

    Bn = 10;
    tPll = 0:Ts:(length(y)*Ts - Ts); % Vector de tiempos completo
    
    % PLL Frecuencia 1s
    [thetaF1] = pll2(y, f1, fs, Bn); % Opcion 1: COSTAS LOOP
    pll_out1 = cos(2*pi*f1*tPll - thetaF1);
    %pll_out1 = pll(y,f1,fs,Bn); % Opcion 2: PLL
    yL_f1 = y.*pll_out1;     
    
    % PLL Frecuencia 0s
    [thetaF0] = pll2(y, f0, fs, Bn); % Opcion 1: COSTAS LOOP
    pll_out0 = sin(2*pi*f0*tPll + thetaF0);
    %pll_out0 = pll(y,f0,fs,Bn); % Opcion 2: PLL
    yL_f0 = y.*pll_out0;

    % Low Pass Filter para Componente de Baja Frecuencia
    N   = 10;       % Order
    Fpass = Rb*1.5;      % Passband Frequency
    Fstop = Rb*4;
    Apass = 0.2;       % Passband Ripple (dB)
    Astop = 15;      % Stopband Attenuation (dB)  

    lpf = designfilt('lowpassfir', 'PassbandFrequency', Fpass, ...
                 'StopbandFrequency', Fstop, 'PassbandRipple', Apass, ...
                 'StopbandAttenuation', Astop, 'SampleRate', fs);
    y0 = filtfilt(lpf,yL_f0);
    y1 = filtfilt(lpf,yL_f1);
    y_dem = (y1 + y0);
    y_dem = 1*y_dem;


end

