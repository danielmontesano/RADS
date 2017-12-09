clear all
close all
clc

%% TX
x = [1 0 1 1 0 1];
% x = ones(1,1000);
Rb = 564.48; %bps

s = moduladorBFSK(x,Rb);

%% CHANNEL

h = 1;
w = 0.01*randn(1,length(s));
y=h.*s+w;

%% RX
[ y_dem ] = demoduladorBFSK(y, Rb, 1);
figure(1)
subplot(3,1,1)
plot(s);
title('Señal Original')
subplot(3,1,2)
plot(y)
title('Señal Recibida')
subplot(3,1,3)
plot(y_dem)
title('Señal Demodulada')

display(['Received Signal: ' num2str(y_dem)])