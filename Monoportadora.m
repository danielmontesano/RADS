clear all
clc

%% TX
x = [0 1 0 1 0 0];
Rb = 564.48; %bps
fs = 56*Rb;

s = moduladorBFSK(x,Rb,fs);

%% CHANNEL
h = 1;
w = 2*randn(1,length(s));
y=h.*s+w;

%% RX
% [ y_dem ] = demoduladorBFSK(y, Rb, 1);
[ y_dem ] = demoduladorBFSK(y, Rb, 3);
% display(['Received Signal: ' num2str(y_dem)])

figure(1)
subplot(3,1,1)
plot(s);
subplot(3,1,2)
plot(y)
subplot(3,1,3)
plot(y_dem)