clear all
close all
clc
sol = [];
%% TX
N = 10000; % Longitud del Mensaje
M = 200; % Numero de simbolos necesarios para el lock del PLL
Z = 200;
% Secuencia de entrenamiento
train = [ones(1,M) zeros(1,M)];
syncS = [];%Start of transmission
for i=1:Z
    if(mod(i,2) == 1)
        syncS = [syncS 1];
    elseif(mod(i,2)==0)
        syncS = [syncS 0];
    end
end
syncE = []; %End of transmission
for i=1:Z/2
    if(mod(i,2) == 1)
        syncE = [syncE 1 1];
    elseif(mod(i,2)==0)
        syncE = [syncE 0 0];
    end
end
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
% h = 1;
% w = 0.02*randn(1,length(s));
% y=h.*s+w;
d = 10000; %Distancia del enlace
disp('Transmitiendo')
y = canalTransmision( s, Rb, fs, fc, d );

%% RX

% Opcion 1: Demodulacion de los simbolos
% Opcion 2: Filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 3: Demodulacion de portadora y filtros paso-banda sintonizados a las frecuencias de los bits con deteccion de energia.
% Opcion 4: PLLs para las frecuencias de los bits.
disp('Demodulando')
[ y_dem ] = demoduladorBFSK(y,Rb,f0,f1,fc,2);
disp('Recuperando')
[ y_sample ] = clockrec( y_dem, 0.1, fs, Rb );

%% RESULTADOS

y_sample(y_sample>=0) = 1;
y_sample(y_sample<0) = 0;

figure
subplot(4,1,1)
plot(s);
title('Sent Signal')
subplot(4,1,2)
plot(y)
title('Received Signal')
subplot(4,1,3)
plot(y_dem)
title('Demoduladated Signal')
subplot(4,1,4)
plot(x);
hold on
plot(y_sample);
title('Original vs Received Signal')

display(['Recovered Signal: ' num2str(y_sample)])

%BER
% [m,i] = max(xcorr(y_sample-0.5,sync-0.5)); %El -0.5 es por xcorr, al sumar 0 no lo hace bien. Sale mejor asi
% start = i-(N+M+M+Z)+Z
start = abs(finddelay(y_sample-0.5,syncS-0.5))+Z+2; %El -0.5 es por xcorr, al sumar 0 no lo hace bien. Sale mejor asi
ends =  abs(finddelay(y_sample-0.5,syncE-0.5));
sizeY = length(y_sample(start:ends))
error = sum((x(M+M+Z+1:end-Z)-y_sample(start:ends)).^2)*100/N; %En porcentaje
figure;
plot(x(M+M+Z+1:M+M+Z+1+20));
hold on
plot(y_sample(start:start+20));

