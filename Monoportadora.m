clear all
close all
clc

%% TX
N = 2000; % Longitud del Mensaje
M = 200; % Numero de simbolos necesarios para el lock del PLL

% Secuencia de entrenamiento
train = [1 0 1 0 ones(1,M) zeros(1,M) 1 0 1];

% Secuencia de entrenamiento + 0s y 1s alternados
x = [ones(1,200) zeros(1,200)];
for i=1:3000
    if(mod(i,2) == 1)
        x = [x 1];
    elseif(mod(i,2)==0)
        x = [x 0];
    end
end

% Secuencia de entrenamiento + 0s y 1s aleatoriamente generados
% x = [ones(1,200) zeros(1,200) round(rand(1,N))];

Rb = 564.48; %bps
fs = 56*Rb;

% Frecuencias de Simbolo
f0 = 8*Rb;
f1 = 7*Rb;

% Frecuancia Portadora y Muestreo
fc = 7.5*Rb;

s = moduladorBFSK(x,Rb,f0,f1,fs);

%% CHANNEL
h = 1;
w = 0.12*randn(1,length(s));
y=h.*s+w;

%% RX

% Opcion 1: Demodulacion de los simbolos
% Opcion 2: Filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 3: Demodulacion de portadora y filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 4: PLLs para las frecuencias de los bits.

[ y_dem ] = demoduladorBFSK(y,Rb,f0,f1,fc,4);
[ y_sample ] = clockrec( y_dem, 0.1, fs, Rb );

%% RESULTADOS
figure
subplot(4,1,1)
plot(x);
title('Original Signal')
subplot(4,1,2)
plot(s);
title('Sent Signal')
subplot(4,1,3)
plot(y)
title('Received Signal')
subplot(4,1,4)
plot(y_dem)
title('Demoduladated Signal')
y_sample(y_sample>=0) = 1;
y_sample(y_sample<0) = 0;

display(['Recovered Signal: ' num2str(y_sample)])
