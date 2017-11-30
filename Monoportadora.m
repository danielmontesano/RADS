clear all
clc

%% TX
x = [1 0 1 1 0 1];
% x = ones(1,1000);
Rb = 564.48; %bps
fs = 56*Rb;

s = moduladorBFSK(x,Rb,fs);

%% CHANNEL

h = 1;
w = 1*randn(1,length(s));
y=h.*s+w;

%% RX
% [ y_dem ] = demoduladorBFSK(y, Rb, 1);
[ y_dem ] = demoduladorBFSK(y, Rb, 4);
% display(['Received Signal: ' num2str(y_dem)])

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