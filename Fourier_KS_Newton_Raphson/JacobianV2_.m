function DF = JacobianV2_(a,n)


% Jacobian matrix for fourier

global L_bar

delta = eye(n+1,n+1);

for k = 1:n+1
    q_k = (k-(n/2)-1)/L_bar;
    for j = 1:n+1

        sum = 0;
        for m = 1:n+1                                            % summation

            if (k-m+n/2+1) <= n && (k-m+n/2+1)>0
                sum = sum + a(m)*delta(k-m+n/2+1,j) + a(k-m+n/2+1)*delta(m,j);
            else
                sum = sum;
            end

        end
        DF(k,j) = (q_k^2 - q_k^4)*delta(k,j) - (1i*q_k/2)*sum;

    end
end

end
