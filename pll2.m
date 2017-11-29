function [theta] = pll2(y, f0, fs, Bn)

fl=3; 
ff=[0 .01 .02 1]; 
fa=[1 1 0 0];
h=firpm(fl,ff,fa);

% Parametros
dseta=1;
T=1/f0;  %sampling frequency

%Computing the P-I loop-filter constants; k1 and k2
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

w0 = 2*pi*f0*T; %frequency (rad/sample)

v_vect = zeros(1, length(y));
theta=zeros(1,length(y)); 

Ts = 1/fs;
Rb = 564.48; %bps
Tb = 1/Rb;
% t = 0:(Ts):(Tb-Ts);
t = Ts:Ts:length(y)*Ts;

% Inicializacion
mu = 0.001;
theta(1)=pi; % estimate vector
zs=zeros(1,fl+1); 
zc=zeros(1,fl+1);
M0=0; M1=0; M2=0;

% PLL
for k=1:length(y)-1
  zs=[zs(2:fl+1), 2*y(k)*sin(2*pi*f0*t(k)+theta(k))];
  zc=[zc(2:fl+1), 2*y(k)*cos(2*pi*f0*t(k)+theta(k))];
  lpfs=fliplr(h)*zs'; 
  lpfc=fliplr(h)*zc'; % output of filters
  theta(k+1)=theta(k) - mu*lpfs*lpfc;   % algorithm update

%   error = lpfs*lpfc;
%   %loop filter
%   v         = error*k1 + error*k2 + M1;
%   M1 = error*k2 + M1;
%   v_vect(k+1) = v;
%   M2 = v + M2 + w0;
%   theta(k+1) = M2;
end

plot(t, theta),
title('Phase Tracking via the Costas Loop')
xlabel('time'); 
ylabel('phase offset')