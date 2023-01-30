clear all;
% Rumayel Hassan Pallock
% This code finds the equilibrium solution of KS equation
% Use sin(0.4x) to get one type of equilibrium
% Use sin(0.8x) to get another type of equilibrium

% Main program
number = 1;
global L
for L = 19.2:0.1:25
    %L = 19.2;
    n = 100;
    x = linspace(-L/2,L/2,n);
    del_x = x(2)-x(1);


    % Initial guess
    if number == 1
        U_0 = sin(.5*x);
        U_0 = U_0';
    else
        U_0 = u_final;
    end
    %U_0 = u_grid(:,450)';  %dedalus data
    

    k = 1;
    err = 10;
    while err> 10e-6
        U_new = U_0 -jacobi(U_0,n)\finite_diff(U_0,n);
        err = max(abs(U_new - U_0));

        U_0 = U_new;
        %plot(x,U_new)
        %pause(0.1)
        if err >10^20
            disp("Solution did not converge")
            break;
        end
    end
    u_final = U_new;
    plot(x,u_final);
    Energy(number) = Calculate_Energy(u_final,del_x,n,L);
    L_bar(number) = L/(2*pi);
    filename = ['L_bar- ' num2str(L_bar(number)) '.jpg'];
    saveas(1,filename);
    number = number + 1
end
save("Energy_Lbar_EQ_2_19_25.mat","Energy","L_bar");

% % Output in a video file
% writerObj = VideoWriter('Video_EQ3.avi');
% writerObj.FrameRate = 10;
% open(writerObj);
% % write the frames to the video
% for i=1:length(f)
%     % convert the image to a frame
%     frame = f(i) ;
%     writeVideo(writerObj, frame);
% end
% close(writerObj);
