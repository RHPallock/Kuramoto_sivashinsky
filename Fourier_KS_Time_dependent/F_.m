function dudt = F_(t,u,kappa,c)
global iter
uhat = fft(u);
duhat = 1i*kappa.*uhat;
dduhat = -(kappa.^2).*uhat;
dddduhat = (kappa.^4).*uhat;

du = ifft(duhat);
ddu = ifft(dduhat);
ddddu = ifft(dddduhat);
dudt = -u.*du - ddu - ddddu;

%iter = iter+1
end
