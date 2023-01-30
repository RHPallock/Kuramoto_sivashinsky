clear all, close all, clc

% Define spatial domain
c = 2;
L = 22;
N = 60;
dx = L/N;
x = -L/2:dx:L/2-dx;

global iter
iter = 1;

% Define discrete wavenumbers
kappa = (2*pi/L)*[-N/2:N/2-1];
kappa = fftshift(kappa');

% Initial condition
u0 = sin(0.45*x);
uhat0 = fft(u0);

%simulate using ode45
dt = 0.01;
t = 0:dt:10000*dt;
[t u] = ode45(@(t,u) F_(t,u,kappa,c),t,u0);

for k = 1:1:length(t)
    plot(x,u(k,:));
    %xlim([-10 10]);
    %ylim([0 1.2]);
    pause(0.1);
    %f(k) = getframe(gcf);
    %filename = ['time- ' num2str(k) '.jpg'];
    %saveas(1,filename);
    
end
% Complex fourier coefficients
%F = fftshift(fft(u(10000,:)))/N;
%m0 = 1+floor(N/2);
%ck = F(m0-N:m0+N);

% Output in a video file
% writerObj = VideoWriter('Video_sin.avi');
% writerObj.FrameRate = 10;
% open(writerObj);
% % write the frames to the video
% for i=1:length(f)
%     % convert the image to a frame
%     frame = f(i) ;    
%     writeVideo(writerObj, frame);
% end
% close(writerObj);


