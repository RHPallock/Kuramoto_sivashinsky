function Func = finite_diff(u,n)
global L

del_x = L/(n-1);

dumm_Func(1) = -u(1)*(u(2)-u(n))/(2*del_x) - (u(2)-2*u(1)+u(n))/(del_x)^2 - (u(3)-4*u(2)+6*u(1)-4*u(n)+u(n-1))/(del_x)^4;
dumm_Func(2) = -u(2)*(u(3)-u(1))/(2*del_x) - (u(3)-2*u(2)+u(1))/(del_x)^2 - (u(4)-4*u(3)+6*u(2)-4*u(1)+u(n))/(del_x)^4;
dumm_Func(n-1) = -u(n-1)*(u(n)-u(n-2))/(2*del_x) - (u(n)-2*u(n-1)+u(n-2))/(del_x)^2 - (u(1)-4*u(n)+6*u(n-1)-4*u(n-2)+u(n-3))/(del_x)^4;
dumm_Func(n) = -u(n)*(u(1)-u(n-1))/(2*del_x) - (u(1)-2*u(n)+u(n-1))/(del_x)^2 - (u(2)-4*u(1)+6*u(n)-4*u(n-1)+u(n-2))/(del_x)^4;
for i = 3:n-2
    dumm_Func(i) = -u(i)*(u(i+1)-u(i-1))/(2*del_x) - (u(i+1)-2*u(i)+u(i-1))/(del_x)^2 - (u(i+2)-4*u(i+1)+6*u(i)-4*u(i-1)+u(i-2))/(del_x)^4;  
end

Func = dumm_Func';
end