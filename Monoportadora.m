clear all
clc

%% TX
x = [0 0 1 0 1 0];
Rb = 564.48; %bps

s = moduladorBFSK(x,Rb);

%% CHANNEL
h = 1;
w = 0.5*randn(1,length(s));
y=h.*s+w;

%% RX
[ y_dem ] = demoduladorBFSK(y, Rb, 1);
display(['Received Signal: ' num2str(y_dem)])