function X_final = Phi(T0,X0)

% if T0 > 0
%     dt = 0.0001;
%     t = 0:dt:T0;
% else
%     dt = 0.0001;
%     t = 0:-dt:T0;
% end

reltol = 1.0e-10; abstol = 1.0e-10;
options = odeset('RelTol',reltol,'AbsTol',abstol);
[T,X] = ode45(@F_,[0,T0],X0,options);

X_final = X(end,:)';

end