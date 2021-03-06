function [cos_out, error_vector] = pll(y, f0, fs, Bn)

% Parametros
kp=1;
k0=1;
dseta=8;
T=1/fs;  %sampling frequency

%Computing the P-I loop-filter constants; k1 and k2
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

%discrete-time simulation 
N         = length(y); %iterations
theta_est = pi/2;
w0        = 2*pi*f0*T;      %frequency (rad/sample)
% deltaw0   = 2*pi*deltaf0*T; %frequency offset (rad/sample)

s_inVector   = zeros(1, N);
error_vector = zeros(1, N);
Real_exp_est = zeros(1, N);
theta_est_vector = zeros(1, N);
cos_out = zeros(1, N);
v_vect = zeros(1, N);
M0=0; M1=0; M2=0;

for n=0:N-1,    
  s_in  =  y(n+1);
  s_est = -sin(theta_est + 2*pi*f0*n*T);
  c_est = cos(theta_est + 2*pi*f0*n*T);
  phase_det = kp*(s_in*s_est);
  error = phase_det;
  
  %plots
  s_inVector(n+1)   = s_in; 
  Real_exp_est(n+1) = s_est;
  error_vector(n+1) = error;
  theta_est_vector(n+1) = theta_est;
  cos_out(n+1) = c_est;
  
  %loop filter
  v         = error*k1 + error*k2 + M1;
  M1 = error*k2 + M1;
  v_vect(n+1) = v;
  M2 = v + M2;
  theta_est_vector(n+1) = v + M2;
  theta_est = M2;
end

% b = 0;
% % 
% subplot(311)
% plot(0:N-1, s_inVector); axis tight; grid
% title('Re\{s_{in}\}')
% 
% figure
% plot(0:N-1, s_inVector, 'b'); hold on
% plot(0:N-1, cos_out, 'r'); axis tight; grid; hold off
% title('Re\{exp_{in}\} (blue) and Re\{s_{est}\} (red)')

% subplot(313)
% plot(0:N-1,v_vect); grid; axis tight

%plot(theta_est_vector)