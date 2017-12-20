clear all
close all
clc

N = 1000; % Longitud del Mensaje
x = round(rand(1,N)); % Mensaje de 0s y 1s de longitud N
Rb = 564.48; %bps
fs = 56*Rb;

order_V = [4.5 5.5 6.5 7.5 8.5 9.5 10.5];
order = order_V(randperm(length(order_V)));
Th = 100; % nº de bits tras los cuales cambia portadora

bloquesFrecuenciales = ceil(N/Th);
y_dem = [];
num = 1;
for i = 1:bloquesFrecuenciales
    fc = order(num)*Rb;
    f0 = (order(num)+0.5)*Rb;
    f1 = (order(num)-0.5)*Rb;
    
    % TX
    start = (i-1)*Th +1;
    fin = i*Th;
    if(fin > length(x))
        fin = length(x);
    end
    s = moduladorBFSK(x(start:fin),Rb,f0,f1,fs);

    % CHANNEL
    h = 1;
    w = 0.32*randn(1,length(s));
    y=h.*s+w;

    % RX
    [ y_dem_i ] = demoduladorBFSK(y,Rb,f0,f1,fc,2);
    y_dem = [y_dem y_dem_i];
    num = num + 1;
    if(num > length(order_V))
        num = num - length(order_V);
    end
end

figure(1)
plot(x,'r')
hold on
plot(y_dem,'b')