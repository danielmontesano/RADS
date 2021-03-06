function [ y, tiempos ] = clockrecExpansion( r, beta, fs, Rb, tiempos, orderIndex )

Ts = 1/fs;
Tb = 1/Rb;
t = 0:(Ts):(Tb-Ts);

l=length(t)/2; 
m=56;
x=r;
% matchfilt=srrc(l,beta,m,0);
% x=conv(r,matchfilt);   
% if(length(r)>2.24e4)
%     tnow=2.24e4;
% else
%     tnow = m*m+1;
% end
tnow = l*2-15;
tau=tiempos(orderIndex); 
tausave(1)=tau; i=0;
mu=2;                            % algorithm stepsize
delta=1;                          % time for derivative
while tnow<length(x)-56              % run iteration
  i=i+1;
  
%   xs(i) = interp1([floor(tnow+tau) (floor(tnow+tau)+1)],...
%                   [x(floor(tnow+tau)) x(floor(tnow+tau)+1)],...
%                   [tnow+tau]);
%   x_deltap = interp1([floor(tnow+tau+delta) (floor(tnow+tau+delta)+1)],...
%                   [x(floor(tnow+tau+delta)) x(floor(tnow+tau+delta)+1)],...
%                   [tnow+tau+delta]);
%   x_deltam = interp1([floor(tnow+tau-delta) (floor(tnow+tau-delta)+1)],...
%                   [x(floor(tnow+tau-delta)) x(floor(tnow+tau-delta)+1)],...
%                   [tnow+tau-delta]);
              
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left

  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-1,1]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

y = xs;

figure
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title(['Constelaci�n de la se�al, f_c = ' num2str(orderIndex+3.5) 'Hz']);
ylabel('S�mbolo estimado' )
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
title(['Offset (\tau), f_c = ' num2str(orderIndex+3.5) 'Hz'])
ylabel('Offset (muestras)'), xlabel('Iteraciones')

tiempos(orderIndex) = tau;

end

