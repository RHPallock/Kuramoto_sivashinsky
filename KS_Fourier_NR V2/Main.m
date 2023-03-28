clear all;
% Newton raphson method for KS equation
% Rumayel Hassan Pallock
% Date: 21/03/2023
format long;
global n L_bar


% Parameters
L = 2*3.1416;
n = 32; % put even numbers of n
dx = L/n;
x = 0:dx:L;
L_bar = 1;

global x_init

% Initial guess
% Load file from TD simulation
a_init = load('TD_0_029910_ID.mat','a');
Init_mod_num = length(a_init.a(:,1));
%X0 = a_init.a(Init_mod_num-10,:)';
X0 = a_init.a(Init_mod_num-10,:)';

min_val = 0.999*ones(n+1,1);
max_val = 1.001*ones(n+1,1);
random_numbers = min_val + rand*(max_val - min_val);
X0 = random_numbers.*X0;
%X0(n/2+1) = [];



%--------------------------
% T guess 
%T0 = .879;
T0 = 0.879;
%--------------------------


flag = 0;
%n = length(X0);
iter = 1;
err(1) = 10;
%STM_perm = STM(T0,X0);
% Algorithm
while err >1e-10
    x_init = X0;

    a = STM_Vectorized(T0,X0) - eye(n+1,n+1);
    a(n/2+1,:) = [];
    a(:,n/2+1) = [];
    b = F_(T0,Phi(T0,X0));
    b(n/2+1,:) = [];
    c_Dumm = F_(T0,X0);
    c_Dumm(n/2+1,:) = [];
    c = transpose(c_Dumm);
    d = 0;

    A = [a b;c d];

    b1_l = -Phi(T0,X0) + X0;
    b1_l(n/2+1,:) = [];
    b2_l = 0;
    B = [b1_l;b2_l];

    Delta = A\B;
    Delta = [Delta(1:n/2);0;Delta(n/2+1:end)];
    %err(iter) = max(abs(Delta(1:end-1,1)));
    err(iter) = max(abs(-Phi(T0,X0) + X0));
    disp("Error = ");
    disp(err(iter))
    X0 = X0 + Delta(1:end-1,1);
    T0 = T0 + Delta(end,1)
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
end
if flag ==0
    disp("Solution found")
end
disp('error')
disp(err)
disp('Time period')
disp(T0)
disp(X0)

% save data into a data file
save('NR_data_trial.mat','X0','T0');

% % Grid points plot
% a_grid = ifftshift(n*X0);
% a_grid = ifft(a_grid);
% plot(x,1i*a_grid);
