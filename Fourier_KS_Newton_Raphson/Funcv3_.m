function dadt = Funcv3_(a,n,L_bar)

%global n L_bar


for k = 1:n+1
    q_k = (k-n/2-1)/L_bar;
    sum = 0;
    for m = 1:n+1
        if  (k-m+n/2+1) <= n+1 && (k-m+n/2+1) >0
            sum = sum + a(m,1)*a(k-m+n/2+1,1);
        else
            sum = sum;

        end
    end
    dadt(k,1) = (q_k^2-q_k^4)*a(k,1) - (1i*q_k/2)*sum;
end

dadt = dadt;
end
