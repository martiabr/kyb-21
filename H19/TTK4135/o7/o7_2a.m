T = 0.1;
k_1 = 1;
k_2 = 1;
k_3 = 1;
A = [1 T; -k_2*T 1-k_1*T];
B = [0; k_3*T];
C = [1 0];
Q = 4*eye(2);
R = 1;
[K,S,e] = dlqr(A,B,0.5*Q,0.5*R)
p = [0.5-0.03i 0.5+0.03i];
L = place(A',C',p).'

n = 50;
N = 60;
t = 0:T:T*n;
x = [5; 1];
x_t = x;
x_hat = [6; 0];
x_hat_t = x_hat;

for i = 1:n
    u_t = -K*x_hat_t;
    x_hat_t = A*x_hat_t + B*u_t + L*C*(x_t - x_hat_t);
    x_t = A*x_t + B*u_t;
    
    x = [x x_t];
    x_hat = [x_hat x_hat_t];
end

phi = [A-B*K B*K; zeros(2) A-L*C]
eig(phi)

plot(t, x_hat(1,:), 'b'),hold,grid;
plot(t, x_hat(2,:), 'b');
plot(t, x(1,:), 'r');
plot(t, x(2,:), 'r');
title('State(r) vs. state estimate(b)');