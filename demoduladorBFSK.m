function [ y_dem ] = demoduladorBFSK(y, Rb, modo)

% Frecuencias de Simbolo
f0 = 8*Rb;
f1 = 7*Rb;

fc=(f0+f1)/2;

% Frecuencia de Muestreo
fs = 56*Rb;
Ts = 1/fs;

Tb = 1/Rb;
t = 0:(Ts):(Tb-Ts);

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

        if(z1>z0)
            a=1;
        elseif(z1<z0)
            a=0;
        end
        y_dem=[y_dem a];
    end
elseif(modo == 2)
    %     F0_l=fc-175-((564.48/2)*(1+0.7));
    %     F0_h=fc+175-((564.48/2)*(1-0.7));
    %     F1_l=fc-175+((564.48/2)*(1-0.7));
    %     F1_h=fc+175+((564.48/2)*(1+0.7));

    F1_l=fc-Rb/2-Rb/2;
    F1_h=fc-Rb/2+Rb/2;
    F0_l=fc+Rb/2-Rb/2;
    F0_h=fc+Rb/2+Rb/2;

    N      = 4;         % Order
    Apass  = 1;         % Passband Ripple (dB)
    Astop  = 10;        % Stopband Attenuation (dB)

    % Construct an FDESIGN object and call its ELLIP method.
    h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, F0_l, F0_h, ...
                          Astop, Apass, Astop, fs);
    bpf_0 = design(h, 'ellip');
    
    h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, F1_l, F1_h, ...
                          Astop, Apass, Astop, fs);
    bpf_1 = design(h, 'ellip');
    
    sLen = length(t);
    y_dem =[];

    for i=sLen:sLen:length(y)
        y0_filt = filter(bpf_0,y((i-sLen+1):i));
        y0=sum(y0_filt.^2);
        y1_filt = filter(bpf_1,y((i-sLen+1):i));
        y1=sum(y1_filt.^2);
%         figure;
%         plot(linspace(0,fs,10^6),10*log10(abs(fft(y0_filt,10^6))));
%         hold on
%         plot(linspace(0,fs,10^6),10*log10(abs(fft(y1_filt,10^6))));
        if(y1<y0)
            a=0;
        elseif(y1>y0)
            a=1;
        end
        y_dem=[y_dem a];
        
     end
    
elseif(modo == 3)   

        
elseif(modo == 4)
    Bn = 25;
    y_dem = [];
    sLen = length(t);    
    [theta] = pll2(y, f0, fs, Bn);
%     for i=sLen:sLen:length(y)
%         [theta] = pll2(y((i-sLen+1):i), f0, fs, Bn);
%         [theta] = pll2(y((i-sLen+1):i), f1, fs, Bn);
%         if(sum(error0.^2)<sum(error1.^2))
%             y_dem = [y_dem 0];
%         else
%             y_dem = [y_dem 1];
%         end
%     end
 end

end


