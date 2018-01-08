function [ y_out ] = canalTransmision( s, Rb, fs, fc, d )
    
    disp_max = 180e-6; %Maxima dispersion
    multi = 1; %Numero de reflexiones de multi trayecto
    
    disp = disp_max*rand(1,multi); %Retardo de las dispersiones
    n_disp = 1000; %ceil(disp*fs); %Array de numero de muestras de dispersion
    
    %Retardo
    n_ret = 1000; %ceil(d*fs/3e8); %numero de muestras en la distancia dada
    y_ret = [zeros(1,n_ret) s zeros(1,max(n_disp))];
    
    
    %Dispersion    
    for k=1:multi
        a = rand(); %Coeficiente de reflexion. Se podria cambiar la distribucion o calcular ademas en funcion del retardo extra
        y_refl = a.*s;
        y_disp = y_ret + [zeros(1,n_ret) zeros(1,n_disp(k)) y_refl zeros(1,max(n_disp)-n_disp(k))];
    end

    
    %Ruido blanco
    w = 0.7 .*randn(1,length(y_disp));

    
    %Interferencias
    Vpk=1;
    Vrms = Vpk/sqrt(2);
    Prms = Vrms^2;
    PSD_s = Prms/Rb;
    sigma = sqrt(PSD_s)/2; %El dos es porque se tiene la misma potencia pero en toda la banda. 
    inter = sigma*randn(1,length(y_disp));
    
    F_l=fc-Rb;
    F_h=fc+Rb;

    N      = 20;      % Order
    Apass  = 1;      % Passband Ripple (dB)
    Astop  = 80;     % Stopband Attenuation (dB)

    
    % Band Pass Filter F1
    h0  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, F_l, F_h, Astop, Apass, Astop, fs);
    bpf = design(h0, 'ellip');
    inter_bp = filter(bpf,inter);
  
    freq = linspace(0,fs, length(y_disp));
    plot(freq, 20*log10(abs(fft(w))));
    hold on;
    plot(freq, 20*log10(abs(fft(y_disp))));
    plot(freq, 20*log10(abs(fft(inter_bp))));
    
    y_out = y_ret + w + inter_bp;

end

