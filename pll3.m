function [theta,fl] = pll3(y, f0, fs, Bn)

fl=500; 
ff=[0 .01 .02 1]; 
fa=[1 1 0 0];
h=firpm(fl,ff,fa);

theta=zeros(1,length(y)); 

Ts = 1/fs;
t = 0:Ts:(length(y)*Ts - Ts);

% Inicializacion
mu = 0.002218;
theta(1)=0.1; % estimate vector
zs=zeros(1,fl+1); 
zc=zeros(1,fl+1);

% PLL
for k=1:length(theta)-1
  zs=[zs(2:fl+1), 2*y(k)*sin(2*pi*f0*t(k)+theta(k))];
  zc=[zc(2:fl+1), 2*y(k)*cos(2*pi*f0*t(k)+theta(k))];
  lpfs=fliplr(h)*zs'; 
  lpfc=fliplr(h)*zc'; % output of filters
  theta(k+1)=theta(k) - mu*lpfs*lpfc;   % algorithm update
end

figure
plot(t, theta),
title('Phase Tracking via the Costas Loop')
xlabel('time'); 
ylabel('phase offset')