function [ y_dem ] = demoduladorBFSK(y, Rb, modo)

% Frecuencias de Simbolo
f0 = 8*Rb;
f1 = 2*Rb;

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
    % A IMPLEMENTAR
elseif(modo == 3)
    % A IMPLEMENTAR   
elseif(modo == 4)
    Bn = 25; 
    y_dem =[];
    sLen = length(t);
    for i=sLen:sLen:length(y)
        [cos_dem0,error0] = pll(y((i-sLen+1):i), f0, fs, Bn);
        [cos_dem1,error1] = pll(y((i-sLen+1):i), f1, fs, Bn);
        if(sum(error0)<sum(error1))
            y_dem = [y_dem cos_dem0];
        else
            y_dem = [y_dem cos_dem1];
        end
    end
end

end

