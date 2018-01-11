clear all
close all
clc

%% TX
N = 10000; % Longitud del Mensaje
M = 0; % Numero de simbolos necesarios para el lock del PLL
Z = 400;
% Secuencia de entrenamiento
train = [ones(1,M) zeros(1,M)];
syncS = round(rand(1,Z));
syncE = round(rand(1,Z));

% Secuencia de entrenamiento + 0s y 1s aleatoriamente generados
x = [train syncS round(rand(1,N)) syncE];

% Secuencia de entrenamiento + 0s y 1s alternados
% x = [ones(1,200) zeros(1,200)];
% for i=1:1000
%     if(mod(i,2) == 1)
%         x = [x 1];
%     elseif(mod(i,2)==0)
%         x = [x 0];
%     end
% end
% x = [x round(rand(1,N))];

Rb = 564.48; %bps
fs = 56*Rb;

% Frecuencias de Simbolo
f0 = 8*Rb;
f1 = 7*Rb;

% Frecuancia Portadora y Muestreo
fc = 7.5*Rb;
disp('Modulando')
s = moduladorBFSK(x,Rb,f0,f1,fs);


%% CHANNEL
d = 10000; %Distancia del enlace
disp('Transmitiendo')
y = canalTransmision( s, Rb, fs, fc, d );
% y = s + 0.01.*randn(1,length(s));

%% RX
% Opcion 1: Filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 2: Demodulacion de portadora y filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 3: PLLs para las frecuencias de los bits.

disp('Demodulando')
[ y_dem ] = demoduladorBFSK(y,Rb,f0,f1,fc,2);
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
subplot(4,1,2)
plot(t, y)
title('Received Signal')
subplot(4,1,3)
plot(t, y_dem)
title('Demoduladated Signal')
subplot(4,1,4)
plot(x);
hold on
plot(y_sample);
title('Original vs Received Signal')

display(['Recovered Signal: ' num2str(y_sample)])

%% BER
nStart = abs(finddelay(y_sample-0.5,syncS-0.5))+Z+1; %El -0.5 es por xcorr, al sumar 0 no lo hace bien. Sale mejor asi
nEnd =  abs(finddelay(y_sample-0.5,syncE-0.5));
%sizeY = length(y_sample(nStart:nEnd));
error = sum((x(M+M+Z+1:end-Z)-y_sample(nStart:nEnd)).^2)*100/N; %En porcentaje
figure;
plot(x(M+M+Z+1:M+M+Z+1+100));
hold on
plot(y_sample(nStart:nStart+100));

% [~,start] = max(xcorr(syncS-0.5,y_sample-0.5));
% [~,ends] = max(xcorr(syncE-0.5,y_sample-0.5));
% nStart = start+2*M;
% nEnds = ends-M;
% lon = length(y_sample);
% error = sum((x - y_sample(lon-nStart+1:end-nEnds)).^2)*100/N  % En porcentaje
% sol = [sol error];
% figure
% plot(y_sample(lon-nStart+1:end-nEnds));
% hold on
% plot(x);
% title('Original vs Received Signal')
