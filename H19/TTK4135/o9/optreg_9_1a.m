rho = 0.5;
c = 0.5;
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
f_grad = @(x) [-400*(x(2)-x(1)^2)*x(1) - 2*(1-x(1)); 200*(x(2)-x(1)^2)];
f_hessian = @(x) [(1200*x(1)^2 - 400*x(2) + 2) -400*x(1); -400*x(1) 200];

H_k = eye(2);
epsilon = 0.001;
x_k = [-1.2; 1];
x_trajectory = [x_k; f(x_k)];
alpha_values = [];


while norm(f_grad(x_k)) > epsilon
    p_k = - f_hessian(x_k)\f_grad(x_k);  %Newton
    %p_k = - f_grad(x_k);  %Steepest
    %p_k = - H_k * f_grad(x_k); %BFGS
    
    alpha = 1;
    while f(x_k + alpha*p_k) > f(x_k) + c*alpha*f_grad(x_k)'*p_k
        alpha = rho*alpha;
    end
    alpha_values = [alpha_values; alpha];
    x_k_new = x_k + alpha*p_k;
    
    s_k = x_k_new - x_k;
    y_k = f_grad(x_k_new) - f_grad(x_k);
    g_k = 1/(y_k'*s_k);
    H_k = (eye(2) - g_k*s_k*y_k')*H_k*(eye(2) - g_k*y_k*s_k') + g_k*s_k*s_k';
    
    x_k = x_k_new;
    x_trajectory = [x_trajectory [x_k; f(x_k)]];
end

disp(x_k)

x = -1.5:0.1:1.5;
y = -1.5:0.1:1.5;
[X,Y] = meshgrid(x,y);
Z = 100*(Y-X.^2).^2+(1-X).^2;
mesh(X,Y,Z);
hold on;
plot3(x_trajectory(1,:), x_trajectory(2,:), x_trajectory(3,:), 'LineWidth', 5);
title('BFGS x0=(-1.2,1)');

figure(2);
plot(x_trajectory(3,:));
grid on;
title('f(k) for BFGS');

figure(3);
plot(alpha_values);
grid on;
title('alpha for each iteration for steepest descent');


