function [ y_dem ] = demoduladorBFSK(y, Rb,f0,f1,fc,modo)

% Frecuencia de Muestreo
fs = 56*Rb;

Ts = 1/fs;

Tb = 1/Rb;
t = 0:(Ts):(Tb-Ts);

sLen = length(t);

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


N      = 20;      % Order
Apass  = 0.1;      % Passband Ripple (dB)
Astop  = 80;     % Stopband Attenuation (dB)

  % Construct an FDESIGN object and call its ELLIP method.
 h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, F0_l, F0_h, Astop, Apass, Astop, fs);
    bpf_0 = design(h, 'ellip');
    
    h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, F1_l, F1_h, ...
                          Astop, Apass, Astop, fs);
    bpf_1 = design(h, 'ellip');
    
    
    y0 = filter(bpf_0,y);
    y1 = filter(bpf_1,y);
    plot(y0)
    
    y0 = y0.^2;
    y1 = y1.^2;
    
    Fpass = Rb;      % Passband Frequency
  

    h = fdesign.lowpass('n,fp,ap,ast', N, Fpass, Apass, Astop, fs);
    bpf = design(h, 'ellip');
    
    y0 = filter(bpf,y0);
    y1 = filter(bpf,y1);
    
    y_dem = y1 - y0;
    
    rect = [zeros(1,sLen) ones(1,sLen) zeros(1,sLen)];
    sync =  srrc(sLen*2, 0.3, 1, 0);
    y_conv = conv(y_dem,sync); %NO MEJORA
    
    
%    y_dem =[];
    
% 
%     for i=sLen:sLen:length(y)
%         y0_filt = filter(bpf_0,y((i-sLen+1):i));
%         y0=sum(y0_filt.^2);
%         y1_filt = filter(bpf_1,y((i-sLen+1):i));
%         y1=sum(y1_filt.^2);
% %         figure;
% %         plot(linspace(0,fs,10^6),10*log10(abs(fft(y0_filt,10^6))));
% %         hold on
% %         plot(linspace(0,fs,10^6),10*log10(abs(fft(y1_filt,10^6))));
%         if(y1<y0)
%             a=0;
%         elseif(y1>y0)
%             a=1;
%         end
%         y_dem=[y_dem a];
%     end

elseif(modo == 3)
    y_dem =[];
    x = Ts:Ts:length(y)/fs;
    LO_1 = cos(2*pi*f1*x);
    LO_0 = cos(2*pi*f0*x);

    dem1 = y.*LO_1;
    dem0 = y.*LO_0;

    N     = 4;       % Order
    Fpass = Rb;      % Passband Frequency
    Apass = 1;       % Passband Ripple (dB)
    Astop = 80;      % Stopband Attenuation (dB)

    h = fdesign.lowpass('n,fp,ap,ast', N, Fpass, Apass, Astop, fs);
    bpf = design(h, 'ellip');
    dem0_filt = filter(bpf,dem0);  
    dem1_filt = filter(bpf,dem1); 
    
    for k=Tb*fs:Tb*fs:length(y)
        if (dem1_filt(k)>dem0_filt(k))
            a=0;
        elseif(dem1_filt(k)<dem0_filt(k))
            a=1;
        end
        y_dem=[y_dem a];       
    end
    
elseif(modo == 4)
    Bn = 50;
    y_dem = [];
    sLen = length(t);  
    [theta] = pll2(y, f1, fs, Bn);
    y_dem = y;
%     for i=sLen:sLen:length(y)
%         [theta,error0] = pll3(y((i-sLen+1):i), f0, fs, Bn);
%         [theta,error1] = pll3(y((i-sLen+1):i), f1, fs, Bn);
%         if(sum(error0.^2)<sum(error1.^2))
%             y_dem = [y_dem 0];
%         else
%             y_dem = [y_dem 1];
%         end
%     end
 end

end


