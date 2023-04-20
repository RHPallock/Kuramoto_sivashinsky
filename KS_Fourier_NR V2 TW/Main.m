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

%global x_init

% Initial guess
% Load file from TD simulation
%a_init = load('equi_2.mat');
%a_init = load("equi_1.mat");
a_init = load('TW_2_initial.mat','a_coeff','a_grid');

%min_val = 1*ones(n+1,1);
%max_val = 1.5*ones(n+1,1);
%random_numbers = min_val + rand*(max_val - min_val);
%X0 = (random_numbers.*a_init.u_init);
X0 = (a_init.a_coeff);
%X0 = a_init.a(1,:)';
%plot(x,a_init.a_grid);
X0(n/2+1) = 0;

%plot(x,a_init.u_final);
plot(x,ifft(ifftshift((n+1)*X0)));
pause(0.1);





%--------------------------
% T guess (Use a fixed T)
%T0 = 1;
T0 =0.5;
%--------------------------
% Velocity guess
%s0 = 0.34;
%TW_2
s0 = -0.12;
% Initial velocity guess is s0*L_bar
disp('Initial velocity guess');
disp(s0*L_bar)

c0_x = s0*T0;


flag = 0;
%n = length(X0);
iter = 1;
err(1) = 10;
%STM_perm = STM(T0,X0);
% Algorithm
while err >1e-10
%     x_init = X0;

    k = linspace(-(n)/2,(n)/2,n+1);
    A_c = diag(exp(1i*k*c0_x));
%     A_c = zeros(n+1,n+1);
%     for j = 1:n+1
%         A_c(j,j) = exp(1i*k(j)*c0);
%     end

    a = STM_Vectorized(T0,X0) - A_c;
    a(n/2+1,:) = [];
    a(:,n/2+1) = [];

    delA_delc = diag((1i*k).*(exp(1i*k*c0_x)));
%     delA_delc = zeros(n+1,n+1);
%     for j = 1:n+1
%         delA_delc(j,j) = (1i*k(j))*exp(1i*k(j)*c0);
%     end

    b = -delA_delc*X0;
    b(n/2+1,:) = [];
    c_Dumm = F_(T0,X0);
    c_Dumm(n/2+1,:) = [];
    c = transpose(c_Dumm);
    d = 0;

    A = [a b;c d];

    b1_l = -(Phi(T0,X0)) + A_c*X0;
    b1_l(n/2+1,:) = [];
    b2_l = 0;
    B = [b1_l;b2_l];

    Delta = A\B;
    Delta = [Delta(1:n/2);0;Delta(n/2+1:end)];
    err(iter) = max(abs(-(Phi(T0,X0)) + A_c*X0));
    disp("Error = ");
    disp(err(iter))
    X0 = X0 + Delta(1:end-1,1);
    c0_x = c0_x + (Delta(end,1));
    (c0_x/T0)*(L_bar)
    %if iter > 20 && err(iter) > 100
    %    X0 = Phi(-T0,X0)
    %end
    if  err(iter) > 10000 
        disp("Solution did not converge")
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
    a_grid = ifftshift((n+1)*X0);
    a_grid = ifft(a_grid);
    plot(x,a_grid)
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
disp((c0_x/T0)*(L_bar))
%disp(X0)

% save data into a data file
save('NR_data_trial.mat','X0','T0');

% % Grid points plot
% a_grid = ifftshift(n*X0);
% a_grid = ifft(a_grid);
% plot(x,1i*a_grid);
