function Jaco = jacobi(u,n)
global L;

del_x = L/(n-1);

A(1,1) = -(u(2)-u(n))/(2*del_x) + 2/del_x^2 - 6/del_x^4;
A(1,n) = u(1)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(1,2) = -u(1)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(1,3) = -1/del_x^4;
A(1,n-1) = -1/del_x^4;

A(2,2) = -(u(3)-u(1))/(2*del_x) + 2/del_x^2 - 6/del_x^4;
A(2,1) = u(2)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(2,3) = -u(2)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(2,4) = -1/del_x^4;
A(2,n) = -1/del_x^4;

A(n,n) = -(u(1)-u(n-1))/(2*del_x) + 2/del_x^2 - 6/del_x^4;
A(n,n-1) = u(n)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(n,1) = -u(n)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(n,2) = -1/del_x^4;
A(n,n-2) = -1/del_x^4;

A(n-1,n-1) = -(u(n)-u(n-2))/(2*del_x) + 2/del_x^2 - 6/del_x^4;
A(n-1,n-2) = u(n-1)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(n-1,n) = -u(n-1)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
A(n-1,1) = -1/del_x^4;
A(n-1,n-3) = -1/del_x^4;

for i = 3:n-2
    for j = 3:n-2
    if i == j
        A(i,j) = -(u(i+1)-u(i-1))/(2*del_x) + 2/del_x^2 - 6/del_x^4;
    elseif j == i-1
        A(i,j) = u(i)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
    elseif j == i+1
        A(i,j) = -u(i)/(2*del_x) - 1/del_x^2 + 4/del_x^4;
    elseif j == i+2
        A(i,j) = -1/del_x^4;
    elseif j == i-2
        A(i,j) = -1/del_x^4;
    else
        A(i,j) = 0;
    end
    end
    
end
Jaco = A;
end

