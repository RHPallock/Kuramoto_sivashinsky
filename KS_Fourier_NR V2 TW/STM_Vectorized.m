function STM_output = STM_Vectorized(T0,X0)
global n
I = eye(n+1,n+1);

% Initial condition with EYE matrix
I = I(:);
I = cat(1,X0,I);


% Integration
reltol = 1.0e-10; abstol = 1.0e-10;
options = odeset('RelTol',reltol,'AbsTol',abstol);
[t,phi] = ode45(@(t,phi) JacobianV2_(t,phi,n),[0,T0],I,options);
nt = (n+1)*(n+1);

% Delete first n+1 entries of phi. Because they are solution of RHS of the
% equations
phi(:,1:n+1) = [];


% Only take the last entries.
STM_output = [phi(end,1:nt/(n+1)).'];
for k = 1:n
    STM_output = cat(2,STM_output,phi(end,((k*nt/(n+1))+1):(k+1)*(nt/(n+1))).');
end
% As ode45 makes imaginary numbers out of phase by pi. So we need to take
% the complex conjugate
%STM_output = conj(STM_output);
end