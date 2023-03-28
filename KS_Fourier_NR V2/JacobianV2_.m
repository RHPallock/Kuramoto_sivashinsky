function DF = JacobianV2_(T0,A0,n)
global Traj_time Guess_traj

%First n+1 number of variables are for the RHS.
x1_ = A0(1:n+1);
x2_ = A0(n+2:end);

% 1st part
% RHS of the dot equation -----------------------------
a = x1_;
for k = 1:n+1
    q_k = (k-n/2-1);
    sum = 0;
    for m = 1:n+1
        if  (k-m+n/2+1) <= n+1 && (k-m+n/2+1) >0
            sum = sum + a(m,1)*a(k-m+n/2+1,1);
        else
            sum = sum;
        end
    end
    %dadt(k,1) = (q_k^2-.029910*q_k^4)*a(k,1) + (1i*q_k)*sum; 
    dadt(k,1) = (q_k^2-0.029910*q_k^4)*a(k,1) - (q_k)*sum; 
end
x_t = dadt(:);

% -----------------------------------------------------------

% 2nd part, use initial x1_
a = x1_;

% Required for vectorization: Converts n^2x1 vector to nxn matrix in this
% case
x2_ = reshape(x2_,[],n+1);

% Defining dirac-delta function
delta = eye(n+1,n+1);

% Matrix formation, calculating Jacobian at each point of guess trajectory.
for k = 1:n+1
    q_k = (k-(n/2)-1);
    for j = 1:n+1

        sum = 0;
        for m = 1:n+1                                            % summation

            if (k-m+n/2+1) <= n && (k-m+n/2+1)>0
                sum = sum + a(m)*delta(k-m+n/2+1,j) + a(k-m+n/2+1)*delta(m,j);
            else
                sum = sum;
            end

        end
        DF(k,j) = (q_k^2 - 0.029910*q_k^4)*delta(k,j) - (q_k)*sum;
        %DF(k,j) = (q_k^2 - 0.029910*q_k^4)*delta(k,j) + 1i*(q_k)*sum;

    end
end
DF = DF*x2_;
DF = DF(:);

% Concatenate solution of first part and 2nd part
DF = cat(1,x_t,DF);

end