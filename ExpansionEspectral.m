clear all
close all
clc

addpath('./Monoportadora')
addpath('./Expansion_Espectral')

N = 10000; % Longitud del Mensaje

%% Secuencia de entrenamiento
x = round(rand(1,N)); % Mensaje de 0s y 1s de longitud N

Rb = 564.48; %bps
Tb = 1/Rb;
fs = 56*Rb;
Ts = 1/fs;

order_V = [4.5 5.5 6.5 7.5 8.5 9.5 10.5];
order = order_V(randperm(length(order_V)));
Th = 500; % n� de bits tras los cuales cambia portadora
t = 0:(Ts):(Tb-Ts);
sLen = length(t);

bloquesFrecuenciales = ceil(N/Th);
s = [];
y_sample = [];
num = 1;

%% TX
fases = zeros(2,length(order_V));
tiempos = zeros(2,length(order_V));
fase0 = 0.8;
for i = 1:bloquesFrecuenciales
    fc = order(num)*Rb;
    f0 = (order(num)+0.5)*Rb;
    f1 = (order(num)-0.5)*Rb;
    
    
    start = (i-1)*Th +1;
    fin = i*Th;
    if(fin > length(x))
        fin = length(x);
    end
    % Se almacena la fase para asegurar continuidad en fase
    [s_i,fase0] = moduladorExpansion(x(start:fin),Rb,f0,f1,fs,fase0);
    s = [s s_i];
    num = num + 1;
    if(num > length(order_V))
        num = num - length(order_V);
    end
end

%% CHANNEL
h = 1;
w = 0.01*randn(1,length(s));
y=h.*s+w;

%% Rx
num = 1;
start = 0; fin = 0;
for i = 1:bloquesFrecuenciales
    fc = order(num)*Rb;
    f0 = (order(num)+0.5)*Rb;
    f1 = (order(num)-0.5)*Rb;
    orderIndex = order(num)-3.5;
    % RX
    start = fin + 1;
    fin = fin + Th*sLen +1;
    [ y_dem_i,fases ] = demoduladorExpansion(y(start:fin),Rb,f0,f1,fases,orderIndex);
    [ y_sample_i, tiempos ] = clockRecExpansion( y_dem_i, 0.1, fs, Rb, tiempos, orderIndex );
    y_sample = [y_sample y_sample_i];
    num = num + 1;
    if(num > length(order_V))
        num = num - length(order_V);
    end
end

y_sample(y_sample>=0) = 1;
y_sample(y_sample<0) = 0;

figure
plot(x,'r')
hold on
plot(y_sample,'b')