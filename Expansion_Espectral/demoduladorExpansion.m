function [ y_dem, fases ] = demoduladorExpansion( y, Rb,f0,f1,fases,orderNum)

    % Frecuencia de Muestreo
    fs = 56*Rb;
    Ts = 1/fs;
    Tb = 1/Rb;

    Bn = 10;
    tPll = 0:Ts:(length(y)*Ts - Ts); % Vector de tiempos completo
    
    % PLL Frecuencia 1s
    [pll_out1,fases(1,orderNum)]  = pllExpansion(y,f1,fs,Bn,fases(1,orderNum)); % Opcion 2: PLL
    yL_f1 = y.*pll_out1;     
    
    % PLL Frecuencia 0s
    [pll_out0,fases(2,orderNum)] = pllExpansion(y,f0,fs,Bn,fases(2,orderNum)); % Opcion 2: PLL
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
    y_dem = (y1 - y0);
    y_dem = 1*y_dem;


end

