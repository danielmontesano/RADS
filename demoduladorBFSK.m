function [ y_dem ] = demoduladorBFSK(y, Rb,f0,f1,fc,modo)

% Frecuencia de Muestreo
fs = 56*Rb;
Ts = 1/fs;
Tb = 1/Rb;

%% MODO 1: FILTROS PASO BANDA SINTONIZADOS  
if(modo == 1)
    F1_l=fc-Rb/2-Rb/2;
    F1_h=fc-Rb/2+Rb/2;
    F0_l=fc+Rb/2-Rb/2;
    F0_h=fc+Rb/2+Rb/2;

    N = 2;      % Order

    % Band Pass Filter F1   
    h0 = fdesign.bandpass('n,f3db1,f3db2', N, F0_l, F0_h, fs);
    bpf_0 = design(h0, 'butter');
       
    y0 = filter(bpf_0,y);
    
    % Band Pass Filter F0
    h1 = fdesign.bandpass('n,f3db1,f3db2', N, F1_l, F1_h, fs);
    bpf_1 = design(h1, 'butter');   
    y1 = filter(bpf_1,y);
       
    % Low Pass Filter para Envolvente
    Fpass = Rb;    
    h = fdesign.lowpass('n,f3db', N, Fpass, fs);
    lpf = design(h, 'butter');
    
    % Obtencion de Envolvente (Simbolos)
    y0 = filter(lpf,y0.^2);
    y1 = filter(lpf,y1.^2);   
    y_dem = (y1 - y0);
    y_dem = 2.5*y_dem;
    
%% MODO 2: DEMODULACIÓN NO COHERENTE Y FILTROS PASO BAJO
elseif(modo == 2)
    x = Ts:Ts:length(y)/fs;
    LO_1 = cos(2*pi*f1*x);
    LO_0 = cos(2*pi*f0*x);
        
    dem1 = y.*LO_1;
    dem0 = y.*LO_0;

    N     = 2;       % Order
    Fpass = Rb;      % Passband Frequency

    h = fdesign.lowpass('n,f3db', N, Fpass, fs);
    lpf = design(h, 'butter');
    dem0_filt = filter(lpf,dem0);  
    dem1_filt = filter(lpf,dem1); 
    
    y0 = filter(lpf,dem0_filt.^2);
    y1 = filter(lpf,dem1_filt.^2);
    y_dem = (y1 - y0);
    y_dem = 6.5*y_dem;

%% MODO 3: UTILIZACION DEL PLL
elseif(modo == 3)
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

end


