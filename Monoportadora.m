clear all
close all
clc

%% TX
N = 100; % Longitud del Mensaje
% x = ones(1,1000); % Solo 1s
x = round(rand(1,N)); % Mensaje de 0s y 1s de longitud N
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
w = 0.2*randn(1,length(s));
y=h.*s+w;

%% RX

% Opcion 1: Demodulaci�n de los simbolos
% Opcion 2: Filtros paso-banda sintonizados a las frecuencias de los bits con detecci�n de energ�a.
% Opcion 3: Demodulaci�n de portadora y filtros paso-banda sintonizados a las frecuencias de los bits con detecci�n de energ�a.
% Opcion 4: PLLs para las frecuencias de los bits.

[ y_dem ] = demoduladorBFSK(y,Rb,f0,f1,fc,2);
figure(1)
subplot(3,1,1)
plot(s);
title('Se�al Original')
subplot(3,1,2)
plot(y)
title('Se�al Recibida')
subplot(3,1,3)
plot(y_dem)
title('Se�al Demodulada')

display(['Received Signal: ' num2str(y_dem)])