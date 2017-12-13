function [ y_out ] = canalTransmision( y, Rb, fs, d )
%CANALTRANSMISION Summary of this function goes here
%   Detailed explanation goes here
    
    disp_max = 180e-6; %Maxima dispersion
    multi = 2; %Numero de reflexiones de multi trayecto
    
    disp = disp_max*rand(1,multi); %Retardo de las dispersiones
    n_disp = ceil(disp*fs); %Array de numero de muestras de dispersion
    
    %Retardo
    n_ret = ceil(d*fs/3e8); %numero de muestras en la distancia dada
    y_ret = [zeros(1,n_ret) y zeros(1,max(n_disp))];
    
    %Dispersion
   
    
    for k=1:multi
        a = rand(); %Coeficiente de reflexion. Se podria cambiar la distribucion o calcular ademas en funcion del retardo extra
        y_refl = a.*y;
        y_disp = y_ret + [zeros(1,n_ret) zeros(1,n_disp(k)) y_refl zeros(1,max(n_disp)-n_disp(k))];
    end
    
    %Ruido blanco
    w = 0.1*randn(1,length(y_ret));
    
    %Interferencias
    
    
    y_out = y_ret + w;
    plot(y_out)
end

