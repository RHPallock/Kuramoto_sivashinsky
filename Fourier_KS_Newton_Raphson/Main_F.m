clear all
% Fourier method
% Rumayel hassan pallock
% Date 1/30/2023


% Parameters
L = 19;
n = 64; % put even numbers of n
dx = L/n;
x = -L/2:dx:L/2;

global iter L_bar
iter = 1;

L_bar = L/(2*pi);



% Initial guess

for j = -n/2:n/2
    c(j+n/2+1) = (1/(L))*integral(@(x) sin(.2*x).*exp(-1i.*pi.*j.*x./L),-L/2,L/2);
end



%plot(x,uhat0);
u_init = c';


err = 20;


while err>10e-6
    F_mod = Funcv3_(u_init,n,L_bar);
    F_mod(n/2+1,:) = [];
    jack_mod = JacobianV2_(u_init,n);
    jack_mod(n/2+1,:) = [];
    jack_mod(:,n/2+1) = [];
    u_init(n/2+1,:) = [];
    u_new = u_init - jack_mod\F_mod;
    err = max(abs(u_new-u_init));
    u_init = u_new;
    u_init = [u_init(1:n/2);0;u_init(n/2+1:end)];
    if err > 10^20
         disp('Solution did not converge');
         break;
    end
    plot(x,real(u_init));
    hold on
    plot(x,imag(u_init));
    hold off;
    pause(0.1);
    legend('Cos coefficient','sine coefficient');
end

% Going to the grid mode
sum = 0;
for j = -n/2:n/2
    sum = sum + 1*u_init(j+n/2+1)*exp(1i*2*pi*j*x'/(L));
end
figure(2);
plot(x,sum)
