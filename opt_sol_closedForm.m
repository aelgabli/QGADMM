function out=opt_sol_closedForm(XX,YY)
%s = size;
%cvx_begin %quiet
%cvx_precision low
%cvx_solver mosek
%variable x(s)
%minimize(0.5*sum_square(XX*x - YY))
x = (XX'*XX)\(XX'*YY);
%cvx_end
out =0.5*norm(XX*x - YY)^2;
% clear
% close all
% load anis.mat
% %load anis2.mat
% 
% cvx_begin
% cvx_precision low
% cvx_solver mosek
% variables z1(50)
% minimize(0.5*sum_square(X_fede*z1 - y_fede))
% cvx_end
% out =z1;
% obj0 = (0.5*sum_square(X_fede*z1 - y_fede))
% 
% cvx_begin
% cvx_precision low
% cvx_solver mosek
% %variables z1(50) z2(50)
% variables z2(50)
% minimize(0.5*sum_square(H1*z2 - y1)+0.5*sum_square(H2*z2 - y2))
% %subject to
% %z1 == z2
% cvx_end
% out =z2;
% obj1 = (0.5*sum_square(H1*z2 - y1)+0.5*sum_square(H2*z2 - y2))
% 
% loss=0.5*norm(X_fede*z2-y_fede,2)^2;
% 
% result = loss-loss2(end)
% 
% lambda = zeros(50,1);
% x1 = zeros(50,1);
% x2 = zeros(50,1);
% rho =0.1;
% max_iter = 50;
%  for i = 1:max_iter
% cvx_begin quiet
% cvx_precision low
% cvx_solver mosek
% variable x(50)
% minimize(0.5*sum_square(H1*x - y1)+ lambda'*x + rho/2*sum_square(x-x2))
% cvx_end
% x1 =x;
% 
% cvx_begin quiet
% cvx_precision low
% cvx_solver mosek
% variable x(50)
% minimize(0.5*sum_square(H2*x - y2)- lambda'*x + rho/2*sum_square(x-x1))
% cvx_end
% x2 =x;
% 
% obj2 = (0.5*sum_square(H1*x1 - y1)+0.5*sum_square(H2*x2 - y2))
% lambda = lambda + rho*(x1-x2);
% 
% loss=0.5*norm(X_fede*x1-y_fede,2)^2;
% 
% result = loss-loss2(end);
%     
% 
% end
     




