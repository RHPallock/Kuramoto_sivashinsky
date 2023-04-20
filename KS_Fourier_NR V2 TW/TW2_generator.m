clear all;
a0 =     0.03171;
a1 =    -0.01544;
b1 =       -0.86;
a2 =      0.3255;
b2 =      0.3691;
a3 =       0.847;
b3 =     -0.4846;
a4 =      0.2365;
b4 =     -0.6154;
a5 =      0.0218;
b5 =    -0.02424;
a6 =      0.0705;
b6 =     0.09153;
a7 =     0.07938;
b7 =      0.0314;
a8 =     0.02109;
b8 =    -0.01453;
w =       1;
n = 32;
L = 22;
x = linspace(-L/2,L/2,n+1)';
% Plot the trial fit
f_x = a0 + a1*cos(2*pi*x*w/22) + b1*sin(2*pi*x*w/22) + a2*cos(2*pi*2*x*w/22) + b2*sin(2*pi*2*x*w/22) + a3*cos(2*pi*3*x*w/22) + b3*sin(2*pi*3*x*w/22) + ...
    a4*cos(2*pi*4*x*w/22) + b4*sin(2*pi*4*x*w/22) + a5*cos(2*pi*5*x*w/22) + b5*sin(2*pi*5*x*w/22) + ...
    a6*cos(2*pi*6*x*w/22) + b6*sin(2*pi*6*x*w/22) + a7*cos(2*pi*7*x*w/22) + b7*sin(2*pi*7*x*w/22) + ...
    a8*cos(2*pi*8*x*w/22) + b8*sin(2*pi*8*x*w/22);
figure(1)
plot(x,f_x);
pause(0.1);

a_coeff = fft(f_x);
a_coeff = (1/(n+1))*fftshift(a_coeff);

figure(2);
a_grid = ifftshift((n+1)*a_coeff);
a_grid = ifft(a_grid);
plot(x,a_grid);

save('TW_2_initial.mat','a_coeff','a_grid','x');