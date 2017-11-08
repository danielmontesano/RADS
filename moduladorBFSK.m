% function [ output_args ] = moduladorBFSK( input_args )

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
w = 0;
y=h.*s+w;

%% RX
