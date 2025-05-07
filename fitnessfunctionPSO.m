function [ f ] = fitnessfunctionPSO(e,Fs,Ft )
velocity=1e1;
Fs=Fs/velocity;
Ft=Ft/velocity;
if Fs>Ft
f=2*abs(Fs+(1*(velocity^2)))/1e3;
else
f=2*abs(Ft+(1*(velocity^2)))/1e3;
end
end

