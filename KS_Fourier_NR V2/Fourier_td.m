clear all;
% Time dependent simulation

global n L_bar


% Parameters
L = 2*3.1416;
n = 32; % put even numbers of n
dx = L/n;
x = 0:dx:L;
L_bar = 1;

% Conversion from grid values to fourier coefficient
% Initial condition
% Valid for function
% for j = -n/2:n/2
%     c(j+n/2+1) = (1/(L))*integral(@(Y) 0.01*Y.*exp(-1i.*2.*pi.*j.*Y./L),0,2*L);
% end
%c = fft(0.01*x);
%c = (1/n)*fftshift(c);

c = load('NR_data_trial_PO_1.mat','X0','T0');

%a_init = imag(c');
a_init = c.X0;
a_init(n/2+1) = 0;

% Time guess
T0 = c.T0;

% Time span
if T0 > 0
    dt = 0.0001;
    t = 0:dt:T0;
else
    dt = 0.0001;
    t = 0:-dt:T0;
end



% Integration
reltol = 1.0e-10; abstol = 1.0e-10;
options = odeset('RelTol',reltol,'AbsTol',abstol);
[t,a] = ode45(@Funcv3_,t,a_init,options);

% for time = 1:length(t)
%     plot(x,a(time,:))
%     pause(0.1)
% end


% % Going to the grid mode
x = 0:0.01*dx:L;
figure(1)
for time = 1:100:length(t)
         sum = 0;
         for j = -n/2:n/2
             sum = sum + 1i*a(time,j+n/2+1)*exp(1i*2*pi*j*x'/(L));
         end
    %a_grid = ifftshift(n*a(time,:));
    %a_grid = ifft(a_grid);
    %plot(x,1i*a_grid);
    plot(x,sum);
    pause(0.1)
    %f(time) = getframe(gcf);
    %     filename = ['time-10x ' num2str(time) '.jpg'];
    %     saveas(1,filename);
end

% % % % Output in a video file
%  writerObj = VideoWriter('niu_0046.avi');
%  writerObj.FrameRate = 10;
%  open(writerObj);
%  % write the frames to the video
%  for i=length(f)-50000:100:length(f)
%      % convert the image to a frame
%      frame = f(i) ;
%      writeVideo(writerObj, frame);
%  end
%  close(writerObj);


figure(2)
% Phase portrait
const = 0;
const2 =1;
plot3(a(const2:end,n/2+2),a(const2:end,n/2+3),a(const2:end,n/2+4));
hold on;
scatter3(a(end-const,n/2+2),a(end-const,n/2+3),a(end-const,n/2+4))
hold off;
xlabel('a1');
ylabel('a2');
zlabel('a3');