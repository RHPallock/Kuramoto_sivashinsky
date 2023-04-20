clear all;
% Newton raphson method for KS equation Travelling wave
% Rumayel Hassan Pallock
% Date: 04/03/2023
format long;
global n L_bar


% Parameters
L = 22;
n = 32; % put even numbers of n
dx = L/n;
x = -L/2:dx:L/2;
L_bar = L/(2*pi);

global x_init

% Initial guess
% Load file from TD simulation
a_init = load('equi_2.mat','a');


X0 = a_init.a;
plot(x,ifft(ifftshift((n+1)*X0)));
pause(0.01);





%--------------------------
% T guess (Use a fixed T)
T0 = 1; 
%--------------------------
% Velocity guess
s0 = .1;
c0 = s0*T0;


flag = 0;
%n = length(X0);
iter = 1;
err(1) = 10;
%STM_perm = STM(T0,X0);
% Algorithm
while err >1e-10
    x_init = X0;

    k = linspace(-(n)/2,(n)/2,n+1);
    A_c = diag(exp(1i*k*c0));

    a = STM_Vectorized(T0,X0) - A_c;
    a(n/2+1,:) = [];
    a(:,n/2+1) = [];

    delA_delc = diag((1i*k).*(exp(1i*k*c0)));
    b = -delA_delc*X0;
    b(n/2+1,:) = [];
    c_Dumm = F_(T0,X0);
    c_Dumm(n/2+1,:) = [];
    c = transpose(c_Dumm);
    d = 0;

    A = [a b;c d];

    b1_l = -Phi(T0,X0) + A_c*X0;
    b1_l(n/2+1,:) = [];
    b2_l = 0;
    B = [b1_l;b2_l];

    Delta = A\B;
    Delta = [Delta(1:n/2);0;Delta(n/2+1:end)];
    %err(iter) = max(abs(Delta(1:end-1,1)));
    err(iter) = max(abs(-Phi(T0,X0) + A_c*X0));
    %err(iter) = max(abs(Delta));
    disp("Error = ");
    disp(err(iter))
    X0 = X0 + Delta(1:end-1,1);
    c0 = c0 + real(Delta(end,1))
    %if iter > 20 && err(iter) > 100
    %    X0 = Phi(-T0,X0)
    %end
    if  err(iter) > 10000 || iter > 150
        disp("50 iterations reached")
        flag = 1;
        break;
    end
    if T0<0
        disp('Period is negative');
        break;
    end
    iter = iter + 1
    %scatter(X0(1,1),X0(2,1))
    %hold on;
    a = ifftshift((n+1)*X0);
    a = ifft(a);
    plot(x,a)
    pause(0.01)
end
if flag ==0
    disp("Solution found")
end
disp('error')
disp(err)
disp('Time period')
disp(T0)
disp('velocity');
disp(c/T0)
disp(X0)

% save data into a data file
save('NR_data_trial.mat','X0','T0');

% % Grid points plot
% a_grid = ifftshift(n*X0);
% a_grid = ifft(a_grid);
% plot(x,1i*a_grid);
