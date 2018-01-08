function [ y_dem ] = demoduladorBFSK(y, Rb,f0,f1,fc,modo)

% Frecuencia de Muestreo
fs = 56*Rb;
Ts = 1/fs;
Tb = 1/Rb;
t = 0:(Ts):(Tb-Ts);
sLen = length(t);

%% MODO 1: DEMODULACION SIMPLE
if(modo == 1)
    y_dem =[];
    sLen = length(t);
    for i=sLen:sLen:length(y)
        s_dem1=cos(2*pi*f1*t);
        s_dem0=cos(2*pi*f0*t);

        y_dem1=s_dem1.*y((i-sLen+1):i);
        y_dem0=s_dem0.*y((i-sLen+1):i);

        z1=sum(y_dem1)/length(y_dem1);
        z0=sum(y_dem0)/length(y_dem0);
    %TODO sobra esto
        if(z1>z0)
            a=1;
        elseif(z1<z0)
            a=0;
        end
        y_dem=[y_dem a];
    end
    
%% MODO 2: FILTROS PASO BANDA SINTONIZADOS  
elseif(modo == 2)
    %     F0_l=fc-175-((564.48/2)*(1+0.7));
    %     F0_h=fc+175-((564.48/2)*(1-0.7));
    %     F1_l=fc-175+((564.48/2)*(1-0.7));
    %     F1_h=fc+175+((564.48/2)*(1+0.7));

    F1_l=fc-Rb/2-Rb/2;
    F1_h=fc-Rb/2+Rb/2;
    F0_l=fc+Rb/2-Rb/2;
    F0_h=fc+Rb/2+Rb/2;

    N      = 2;      % Order
    Apass  = 1;      % Passband Ripple (dB)
    Astop  = 80;     % Stopband Attenuation (dB)

    % Band Pass Filter F1
%     h0  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, F0_l, F0_h, Astop, Apass, Astop, fs);
%     bpf_0 = design(h0, 'ellip');
    
%     h0 = fdesign.bandpass('n,f3db1,f3db2', N, F0_l, F0_h, fs);
%     bpf_0 = design(h0, 'butter');
%    bpFilt = designfilt('bandpassfir', ...
%         'CutoffFrequency1',500,'CutoffFrequency2',560, ...
%         'SampleRate',1500);
bpf1 = designfilt('bandpassfir', 'StopbandFrequency1', F1_l*0.8, ...
                    'PassbandFrequency1', F1_l, 'PassbandFrequency2', ...
                    F1_h, 'StopbandFrequency2', F1_l*1.2, ...
                    'StopbandAttenuation1', Astop, 'PassbandRipple', ...
                    Apass, 'StopbandAttenuation2', Astop, 'SampleRate', ...
                    fs);
    y1b = filtfilt(bpf1,y);
    
    bpf0 = designfilt('bandpassfir', 'StopbandFrequency1', F0_l*0.8, ...
                    'PassbandFrequency1', F0_l, 'PassbandFrequency2', ...
                    F0_h, 'StopbandFrequency2', F0_l*1.2, ...
                    'StopbandAttenuation1', Astop, 'PassbandRipple', ...
                    Apass, 'StopbandAttenuation2', Astop, 'SampleRate', ...
                    fs);
    y0b = filtfilt(bpf0,y);
       
       
%     y0 = filter(bpf_0,y);
    
    % Band Pass Filter F0
%     h1  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, F1_l, F1_h, Astop, Apass, Astop, fs);
%     bpf_1 = design(h1, 'ellip');   
%     h1 = fdesign.bandpass('n,f3db1,f3db2', N, F1_l, F1_h, fs);
%     bpf_1 = design(h1, 'butter');
%     
%     y1 = filter(bpf_1,y);
%        
    % Low Pass Filter para Envolvente
%     Fpass = Rb;    
%     h  = fdesign.lowpass('N,Fp,Ap,Ast', N, Fpass, Apass, Astop, fs);
%     lpf = design(h, 'ellip');
%     h = fdesign.lowpass('n,f3db', N, Fpass, fs);
%     lpf = design(h, 'butter');


    Fpass = Rb*2;      % Passband Frequency
    Fstop = Rb*3;
    Apass = 0.1;       % Passband Ripple (dB)
    Astop = 40;      % Stopband Attenuation (dB)

lpf = designfilt('lowpassfir', 'PassbandFrequency', Fpass, ...
                 'StopbandFrequency', Fstop, 'PassbandRipple', Apass, ...
                 'StopbandAttenuation', Astop, 'SampleRate', fs);
     y0 = filtfilt(lpf,y0b.^2);
    y1 = filtfilt(lpf,y1b.^2);
    
    % Obtencion de Envolvente (Simbolos)
%     y0 = filter(lpf,y0.^2);
%     y1 = filter(lpf,y1.^2);   
    y_dem = (y1 - y0);
    y_dem = 2.5*y_dem;
    
%% MODO 3: DEMODULACIÓN NO COHERENTE Y FILTROS PASO BAJO
elseif(modo == 3)
    x = Ts:Ts:length(y)/fs;
    LO_1 = cos(2*pi*f1*x);
    LO_0 = cos(2*pi*f0*x);
    
%     N = 2;
%     Fl = fc-Rb;
%     Fh = fc+Rb;
%     h0 = fdesign.bandpass('n,f3db1,f3db2', N, Fl, Fh, fs);
%     bpf = design(h0, 'butter');
%     y = filter(bpf ,y);
        
    dem1 = y.*LO_1;
    dem0 = y.*LO_0;

    N     = 2;       % Order
    Fpass = Rb;      % Passband Frequency
    Apass = 1;       % Passband Ripple (dB)
    Astop = 80;      % Stopband Attenuation (dB)

%     h = fdesign.lowpass('n,fp,ap,ast', N, Fpass, Apass, Astop, fs);
%     lpf = design(h, 'ellip');
    h = fdesign.lowpass('n,f3db', N, Fpass, fs);
    lpf = design(h, 'butter');
    dem0_filt = filter(lpf,dem0);  
    dem1_filt = filter(lpf,dem1); 
    
    y0 = filter(lpf,dem0_filt.^2);
    y1 = filter(lpf,dem1_filt.^2);
    y_dem = (y1 - y0);
    y_dem = 6.5*y_dem;
    
%     for k=Tb*fs:Tb*fs:length(y)
%         if (dem1_filt(k)>dem0_filt(k))
%             a=0;
%         elseif(dem1_filt(k)<dem0_filt(k))
%             a=1;
%         end
%         y_dem=[y_dem a];       
%     end

%% MODO 4: UTILIZACION DEL PLL
elseif(modo == 4)
    Bn = 10;
    y_dem = [];
    sLen = length(t); % Longitud de cada simbolo    
    tPll = 0:Ts:(length(y)*Ts - Ts); % Vector de tiempos completo
    
    % PLL Frecuencia 1s
    [thetaF1] = pll2(y, f1, fs, Bn); % Opcion 1: COSTAS LOOP
    pll_out1 = cos(2*pi*f1*tPll - thetaF1);
%     pll_out1 = pll(y,f1,fs,Bn); % Opcion 2: PLL
    yL_f1 = y.*pll_out1;     
    
    % PLL Frecuencia 0s
    [thetaF0] = pll2(y, f0, fs, Bn); % Opcion 1: COSTAS LOOP
    pll_out0 = sin(2*pi*f0*tPll + thetaF0);
%     pll_out0 = pll(y,f0,fs,Bn); % Opcion 2: PLL
    yL_f0 = y.*pll_out0;

    % Low Pass Filter para Envolvente
    N     = 4;       % Order
    Fpass = Rb;      % Passband Frequency
    Apass = 1;       % Passband Ripple (dB)
    Astop = 80;      % Stopband Attenuation (dB)
    
     N   = 10;       % Order
    Fpass = Rb*2;      % Passband Frequency
    Fstop = Rb*3;
    Apass = 0.1;       % Passband Ripple (dB)
    Astop = 40;      % Stopband Attenuation (dB)
    match = 'stopband';  % Band to match exactly

%     h  = fdesign.lowpass('N,Fp,Ap,Ast', N, Fpass, Apass, Astop, fs);
%     lpf = design(h, 'ellip');

%     h = fdesign.lowpass('n,f3db', N, Fpass, fs);
%     lpf = design(h, 'butter');
% h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, fs);
% lpf = design(h, 'butter', 'MatchExactly', match);
    
    %TODO parece que no hace falta.
%     yL_f0 = [yL_f0 zeros(1,1000)];
%     yL_f1 = [yL_f1 zeros(1,1000)];
    
    %y0 = filter(lpf,yL_f0);
    %y1 = filter(lpf,yL_f1);

lpf = designfilt('lowpassfir', 'PassbandFrequency', Fpass, ...
                 'StopbandFrequency', Fstop, 'PassbandRipple', Apass, ...
                 'StopbandAttenuation', Astop, 'SampleRate', fs);
     y0 = filtfilt(lpf,yL_f0);
    y1 = filtfilt(lpf,yL_f1);
    y_dem = (y1 + y0);
    y_dem = 1*y_dem;
end

end


