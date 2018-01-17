function [s] = moduladorBFSK(x,Rb,f0,f1,fs)

% Periodos
T0 = 1/f0;
T1 = 1/f1;

% Periodo muestreo
Ts = 1/fs;

Tb = 1/Rb;
t = 0:(Ts):(Tb-Ts);

dPhi1 = (2*pi)/(fs*T1);
dPhi0 = (2*pi)/(fs*T0);

phaseArray = 0.8;

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

end
