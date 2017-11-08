% function [ output_args ] = moduladorBFSK( input_args )

%% TX
x = [0 0 0 1 0 0 1 0 1];
Rb = 563.48;
f0 = 8*Rb;
f1 = 7*Rb;
fc = 7.5*Rb;
fs = 56*Rb;

Tb = 1/Rb;
t = 0:Tb/fs:Tb;
s = [];
phi = 0;
figure
for i = 1:length(x)
    if(x(i) == 1)
        s = [s cos(2*pi*f1*t+ phi)];
    elseif(x(i) == 0)
        s = [s cos(2*pi*f0*t+ phi)];
    end
    plot(s)
    hold on
end

%% CHANNEL
h = 1;
w = 0;
y=h.*s+w;
