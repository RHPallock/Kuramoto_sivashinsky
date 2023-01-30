clear all;
% Rumayel Hassan Pallock
% This code finds the equilibrium solution of KS equation
% Use sin(0.4x) to get one type of equilibrium
% Use sin(0.8x) to get another type of equilibrium

% Main program
number = 1;
global L

L = 22;
n = 100;
x = linspace(-L/2,L/2,n);
del_x = x(2)-x(1);


%Initial guess
U_0 = sin(.3*x);
U_0 = U_0';



k = 1;
err = 10;
while err> 10e-6
    U_new = U_0 -jacobi(U_0,n)\finite_diff(U_0,n);
    err = max(abs(U_new - U_0));

    U_0 = U_new;
    plot(x,U_new)
    pause(0.1)
    if err >10^20
        disp("Solution did not converge")
        break;
    end
end
u_final = U_new;
plot(x,u_final);


