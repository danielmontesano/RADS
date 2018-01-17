clear all
close all
clc
%%Parametros
%Modo de recepcion:
% Opcion 1: Filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 2: Demodulacion de portadora y filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 3: PLLs para las frecuencias de los bits.
modo = 3;

N = 10000; % Longitud del Mensaje
M = 0; % Numero de simbolos necesarios para el lock del PLL
Z = 400; % Numero de simbolos de la secuencia de entrenamiento

d = 10000; %Distancia del enlace
sigma = 1; %Potencia de ruido

%Frecuencia de simbolo
Rb = 564.48; %bps
fs = 56*Rb; %Hz

% Frecuencias de portadoras de cada Simbolo
f0 = 8*Rb; %Hz
f1 = 7*Rb; %Hz

%Frecuencia central
fc = 7.5*Rb; %Hz

%% TX
train = [ones(1,M) zeros(1,M)];
syncS = round(rand(1,Z));
syncE = round(rand(1,Z));

% Secuencia de entrenamiento + 0s y 1s aleatoriamente generados
x = [train syncS round(rand(1,N)) syncE];

% Frecuancia Portadora y Muestreo
disp('Modulando')
s = moduladorBFSK(x,Rb,f0,f1,fs);

%% CHANNEL
disp('Transmitiendo')
y = canalTransmision( s, Rb, fs, fc, d, sigma);

%% RX
disp('Demodulando')
[ y_dem ] = demoduladorBFSK(y,Rb,f0,f1,fc,modo);
disp('Recuperando Reloj')
% figure;
% plot(y_dem)
% hold on;
% plot(detrend(y_dem));
%y_dem = detrend(y_dem);
% y_dem = y_dem./max(abs(y_dem));
[ y_sample ] = clockrec( y_dem, 0.1, fs, Rb );

%% RESULTADOS
y_sample(y_sample>=0) = 1;
y_sample(y_sample<0) = 0;

Ts = 1/fs;

t = 0:Ts:(length(y_dem)*Ts - Ts);

figure
subplot(4,1,1)
plot(s);
title('Sent Signal')
xlabel('Simbolos')
subplot(4,1,2)
plot(t, y)
title('Received Signal')
xlabel('Muestras')
subplot(4,1,3)
plot(t, y_dem)
title('Demoduladated Signal')
xlabel('Muestras')
subplot(4,1,4)
plot(x);
hold on
plot(y_sample);
title('Original vs Received Signal')
xlabel('Simbolos')

display(['Recovered Signal: ' num2str(y_sample)])

%% BER
nStart = abs(finddelay(y_sample-0.5,syncS-0.5))+Z+1; %El -0.5 es por xcorr, al sumar 0 no lo hace bien. Sale mejor asi
nEnd =  abs(finddelay(y_sample-0.5,syncE-0.5));
%sizeY = length(y_sample(nStart:nEnd));
error = sum((x(M+M+Z+1:end-Z)-y_sample(nStart:nEnd)).^2)*100/N; %En porcentaje
fprintf('BER: %.2f %%\n',error)
figure;
plot(x(M+M+Z+1:M+M+Z+1+100));
hold on
plot(y_sample(nStart:nStart+100));
title('First 100 samples of original vs received Signal')
xlabel('Simbolos')
legend('Original signal', 'Received signal')