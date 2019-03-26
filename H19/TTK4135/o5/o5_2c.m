N = 30;
M = 120;
t = (0:M-1);
Q = blkdiag(0,0,2);
R = 2;
G = [kron(eye(N), Q) zeros(N * 3, N); zeros(N, N * 3)  R * eye(N)];
A_d = [0 0 0; 0 0 1; 0.1 -0.79 1.78];
A_real = [0 0 0; 0 0 1; 0.1 -0.855 1.85];
b_d = [1 0 0.1]';
b_real = [1 0 0]';
x_0 = [0 0 1]';
A_eq = [eye(3 * N) + [zeros(3, 3 * N); kron(eye(N - 1), -A_d) zeros(3 * (N - 1), 3)]  kron(eye(N), -b_d)];
b_eq = [A_d * x_0; zeros(3 * (N - 1), 1)];
A = [zeros(2 * N, 3 * N) kron(eye(N), [1; -1])];
b = ones(2 * N, 1);

x1 = [x_0(1)];
x2 = [x_0(2)];
x3 = [x_0(3)];
u = [];
x_k = x_0;
u_k = [];

for i = 1:M
    b_eq = [A_d * x_k; zeros(3 * (N - 1), 1)];
    [z, fval, flag, out] = quadprog(G, [], A, b, A_eq, b_eq);
    u_k = z(3*N + 1);
    u = [u; u_k];
    x_k = A_real * [x1(i); x2(i); x3(i)] + b_real * u_k + normrnd(0,0.2);
    x1 = [x1; x_k(1)];
    x2 = [x2; x_k(2)];
    x3 = [x3; x_k(3)];
end

x3 = x3(1:M);

figure
subplot(211);
stem(t, x3);
grid on
title('y=x3 with r=2 and constraint on u - MPC style - with random noise');
subplot(212);
stairs(t, u);
grid on
title('u with r=2 and constraint on u - MPC style - with random noise');