clear all
close all
clc

%% TX
N = 2000; % Longitud del Mensaje
M = 200; % Numero de simbolos necesarios para el lock del PLL

train = [1 0 1 0 ones(1,M) zeros(1,M) 1 0 1];
x = [train round(rand(1,N))]; % Solo 1s
% x = [ones(1,200) zeros(1,200) 1 0 1 0 1 0 1 1 0 1 0 0 1 0 1 0 1 0 1 0];
x = [ones(1,200) zeros(1,200) round(rand(1,N))];
% x = round(rand(1,N)); % Mensaje de 0s y 1s de longitud N

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
w = 0.02*randn(1,length(s));
y=h.*s+w;

%% RX

% Opcion 1: Demodulaci�n de los simbolos
% Opcion 2: Filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 3: Demodulacion de portadora y filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 4: PLLs para las frecuencias de los bits.

[ y_dem ] = demoduladorBFSK(y,Rb,f0,f1,fc,2);
[ y_clock ] = clockrec( y_dem, 0.1, fs, Rb );

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

display(['Recovered Signal: ' num2str(y_dem)])