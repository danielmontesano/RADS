function [ y_dem ] = demoduladorBFSK(y, Rb, modo)

% Frecuencias de Simbolo
f0 = 8*Rb;
f1 = 2*Rb;

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
%        figure;
%        plot(linspace(0,fs,10^6),10*log10(abs(fft(y0_filt,10^6))));
%        hold on
%        plot(linspace(0,fs,10^6),10*log10(abs(fft(y1_filt,10^6))));

        if(y1<y0)
            a=0
        elseif(y1>y0)
            a=1
        end
        y_dem=[y_dem a];
 end

elseif(modo == 3)
   Bn = 25;
   [LO,error0] = pll(y, fc, fs, Bn);
    figure(2);
     plot(LO);
    dem = LO.*y;
    figure(3);
    plot(dem);

N     = 4;         % Order
Fpass = 1000;      % Passband Frequency
Apass = 1;         % Passband Ripple (dB)
Astop = 80;        % Stopband Attenuation (dB)

h = fdesign.lowpass('n,fp,ap,ast', N, Fpass, Apass, Astop, fs);

bpf = design(h, 'ellip', ...
    'SystemObject', true);

    dem_filt = filter(bpf,dem(1:56));
    figure(4);
    plot(dem_filt);


    
y_dem = 0;

elseif(modo == 4)
    error0acum = [];
    error1acum = [];
    Bn = 25; 
    y_dem =[];
    sLen = length(t);
    for i=sLen:sLen:length(y)
        [cos_dem0,error0] = pll(y((i-sLen+1):i), f0, fs, Bn);
        [cos_dem1,error1] = pll(y((i-sLen+1):i), f1, fs, Bn);
        if(sum(error0.^2)<sum(error1.^2))
            y_dem = [y_dem cos_dem0];
        else
            y_dem = [y_dem cos_dem1];
        end
        error0acum = [error0acum error0];
        error1acum = [error1acum error1];
    end
% figure(2);
% plot(error0acum.^2); hold on; plot(error1acum.^2); hold off;
end






end

