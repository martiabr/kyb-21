%%MPC with input blocking
N = 30;
M = 30;
c_steps = [1,1,2,4,8,14];
N_u = length(c_steps);
t = (0:N-1);
Q = blkdiag(0,0,1);
R = 1;
ulb = -1;
uub = 1;
G = 0.5*[kron(eye(N), Q) zeros(N * 3, N_u); zeros(N_u, N * 3)  R * eye(N_u)];
A_d = [0 0 0; 0 0 1; 0.1 -0.79 1.78];
b_d = [1 0 0.1]';
x_0 = [0 0 1]';

A_eq_l = [eye(3 * N) + [zeros(3, 3 * N); kron(eye(N - 1), -A_d) zeros(3 * (N - 1), 3)]];
A_eq_r = [];
for i = 1:N_u
    A_eq_r = [A_eq_r; zeros(c_steps(i) * 3, i-1), kron(ones(c_steps(i), 1), -b_d), zeros(c_steps(i) * 3, N_u - i)];
end
A_eq = [A_eq_l A_eq_r];
b_eq = [A_d * x_0; zeros(3 * (N - 1), 1)];

A = [zeros(2 * N_u, 3 * N) kron(eye(N_u), [uub; ulb])];
b = ones(2 * N_u, 1);

x1 = [x_0(1)];
x2 = [x_0(2)];
x3 = [x_0(3)];
u = [];
x_k = x_0;
u_k = [];

for i = 1:M
    b_eq = [A_d * x_k; zeros(3 * (N - 1), 1)];
    [z, fval, flag, out] = quadprog(G, [], A, b, A_eq, b_eq);
    x_k = [z(1); z(2); z(3)];
    x1 = [x1; x_k(1)];
    x2 = [x2; x_k(2)];
    x3 = [x3; x_k(3)];
    u = [u; z(3*N + 1)];
end

x3 = x3(1:M);

figure
subplot(211);
stem(t, x3);
grid on
title('y=x3 with r=1');
subplot(212);
stairs(t, u);
grid on
title('u with r=1 and constraint on u with increasing input blocking - MPC style');