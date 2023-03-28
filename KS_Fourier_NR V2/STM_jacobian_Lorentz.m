function x_t = STM_jacobian_Lorentz(T0,X0,n)
global x_init Traj_time Guess_traj



x_guess = interp1(Traj_time,Guess_traj,T0);

% Required for vectorization: Converts 9x1 vector to 3x3 matrix in this
% case
X0 = reshape(X0,[],n);

x = x_guess(1,1);
y = x_guess(1,2);
z = x_guess(1,3);

sigma = 10;
b = 8/3;
r = 28;

A = [-sigma sigma 0;r-z -1 -x;y x -b];

% Vectorized format
x_t = A*X0;
x_t = x_t(:);

% With for loop initial conditions
%x_t = A*X0;
end