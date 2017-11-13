% function [ output_args ] = moduladorBFSK( input_args )
close all
clear all
clc

%% TX
x = [0 0 1 0 1 0];
Rb = 564.48; %bps

% Frecuencias de Simbolo
f0 = 8*Rb;
f1 = 7*Rb;

% Periodos
T0 = 1/f0;
T1 = 1/f1;

% Frecuancia Portadora y Muestreo
fc = 7.5*Rb;
fs = 56*Rb;
Ts = 1/fs;

Tb = 1/Rb;
t = 0:(Ts):(Tb-Ts);
s = [];

dPhi1 = (2*pi)/(fs*T1);
dPhi0 = (2*pi)/(fs*T0);

phaseArray = 0;

for i = 1:length(x)
    if(x(i) == 1)
        p = phaseArray(end) + dPhi1*(1:length(t));
        phaseArray = [phaseArray p];
    elseif(x(i) == 0)
        p = phaseArray(end) + dPhi0*(1:length(t));
        phaseArray = [phaseArray p];
    end
end
phaseArray=wrapTo2Pi(phaseArray);

s = cos(phaseArray);
figure
plot(s)

%% CHANNEL
h = 1;
w = 0; %0.1*randn(1,length(s));
y=h.*s+w;

%% RX
y_dem =[];
sLen = length(t);
for i=sLen:sLen:length(y)
    s_dem1=cos(2*pi*f1*t);
    s_dem0=cos(2*pi*f0*t);
    
    y_dem1=s_dem1.*y((i-sLen+1):i);
    y_dem0=s_dem0.*y((i-sLen+1):i);
    
    z1=trapz(t,y_dem1); 
    z0=trapz(t,y_dem0);
    
    A_dem1=round(2*z1/Tb);
    A_dem0= round(2*z0/Tb);
    
    if(A_dem1>1/2)
        a=1;
    elseif(A_dem0>1/2)
        a=0;
    end
    y_dem=[y_dem a];
end
display(['Received Signal: ' num2str(y_dem)])