function [theta, errorVector] = pll2(y, f0, fs, Bn)

fl=500; 
ff=[0 .01 .02 1]; 
fa=[1 1 0 0];
h=firpm(fl,ff,fa);

% % Parametros
dseta=50;
T=1/fs;  %sampling frequency

%Computing the P-I loop-filter constants; k1 and k2
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

theta=zeros(1,length(y)); 
errorVector=zeros(1,length(y)); 

Ts = 1/fs;
% Rb = 564.48; %bps
% Tb = 1/Rb;
% t = 0:(Ts):(Tb-Ts);
t = 0:Ts:(length(y)*Ts - Ts);

% Inicializacion
theta(1)=0.8; % estimate vector
zs=zeros(1,fl+1); 
zc=zeros(1,fl+1);
M1=0; M2=0;

% PLL
for k=1:length(theta)-1
  zs=[zs(2:fl+1), 2*y(k)*sin(2*pi*f0*t(k)+theta(k))];
  zc=[zc(2:fl+1), 2*y(k)*cos(2*pi*f0*t(k)+theta(k))];
  lpfs=fliplr(h)*zs'; 
  lpfc=fliplr(h)*zc'; % output of filters

  error = lpfs*lpfc;
  v         = error*k1 + error*k2 + M1;
  M1 = error*k2 + M1;
  M2 = v + M2;
  theta(k+1) = M2;
  errorVector(k+1) = error;
end

figure
plot(t, theta),
title('Phase Tracking via the Costas Loop')
xlabel('time'); 
ylabel('phase offset')