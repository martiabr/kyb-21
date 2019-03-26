close all;
T = 0.1;
k_1 = 1;
k_2 = 1;
k_3 = 1;
A = [1 T; -k_2*T 1-k_1*T];
B = [0; k_3*T];
Q = 2*eye(2);
R = 1;

n = 50;
N = 10;
t = 0:T:T*(n-1);
x = [5; 1];
x_t = x;
u_t = 0;
u = [];

G = [kron(eye(N), Q) zeros(N * 2, N); zeros(N, N * 2)  R * eye(N)];
A_eq = [eye(2 * N) + [zeros(2, 2 * N); kron(eye(N - 1), -A) zeros(2 * (N - 1), 2)]  kron(eye(N), -B)];
b_eq = [A * x; zeros(2 * (N - 1), 1)];
A_ineq = [zeros(2 * N, 2 * N) kron(eye(N), [1; -1])];
b_ineq = 4 * ones(2 * N, 1);

for i = 1:n
    b_eq = [A * x_t; zeros(2 * (N - 1), 1)];  % use current estimate as initial value in MPC
    [z, fval, flag, out] = quadprog(G, [], A_ineq, b_ineq, A_eq, b_eq);  % Solve optimalization problem
    
    u_t = z(2*N + 1);  % Get first control input
    
    x_t = A*x_t + B*u_t;  % Apply control input to system
    
    u = [u u_t];
    x = [x x_t];
end

x = x(1:2, 1:n);
figure(1);
plot(t, x(1,:), 'r'), hold, grid;
plot(t, x(2,:), 'r');
title('State without observer');
figure(2)
plot(t, u), grid;