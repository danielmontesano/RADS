function [ y ] = clockrec( r, beta, fs, Rb )

Ts = 1/fs;
Tb = 1/Rb;
t = 0:(Ts):(Tb-Ts);

l=length(t); 
m=56;
x=r;
% matchfilt=srrc(l,beta,m,0);
% x=conv(r,matchfilt);   
% if(length(r)>2.24e4)
%     tnow=2.24e4;
% else
%     tnow = m*m+1;
% end
tnow = round(m/2);
tau=0; 
tausave(1)=tau; i=0;
mu=1.5;                            % algorithm stepsize
delta=1;                          % time for derivative
while tnow<length(x)              % run iteration
  i=i+1;
  
  xs(i) = interp1([floor(tnow+tau) (floor(tnow+tau)+1)],...
                  [x(floor(tnow+tau)) x(floor(tnow+tau)+1)],...
                  [tnow+tau]);
  x_deltap = interp1([floor(tnow+tau+delta) (floor(tnow+tau+delta)+1)],...
                  [x(floor(tnow+tau+delta)) x(floor(tnow+tau+delta)+1)],...
                  [tnow+tau+delta]);
  x_deltam = interp1([floor(tnow+tau-delta) (floor(tnow+tau-delta)+1)],...
                  [x(floor(tnow+tau-delta)) x(floor(tnow+tau-delta)+1)],...
                  [tnow+tau-delta]);
              
%   xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
%   x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
%   x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left

  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-1,1]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

y = xs;

figure
subplot(3,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(3,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')
subplot(3,1,3)
plot(r)

end

