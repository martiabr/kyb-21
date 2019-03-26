N = 30;
Q = blkdiag(0,0,2);
R = 15;
G = [kron(eye(N), Q) zeros(N * 3, N); zeros(N, N * 3)  R * eye(N)];
A = [0 0 0; 0 0 1; 0.1 -0.79 1.78];
b = [1 0 0.1]';
x_0 = [0 0 1]';
A_eq = [eye(3 * N) + [zeros(3, 3 * N); kron(eye(N - 1), -A) zeros(3 * (N - 1), 3)]  kron(eye(N), -b)];
b_eq = [A * x_0; zeros(3 * (N - 1), 1)];
%KKT_matrix = [G -A_eq'; A_eq zeros(3 * N)];
%w = KKT_matrix \ [zeros(4 * N, 1); b_eq];
[w, fval, flag, out] = quadprog(G, [], [], [], A_eq, b_eq);
out.iterations
x = w(1:3*N);
u = w(3*N+1:4*N);
x1 = [x_0(1)];
x2 = [x_0(2)];
x3 = [x_0(3)];

for i = 1:N
    x1 = [x1; x(3*(i-1) + 1)];
    x2 = [x2; x(3*(i-1) + 2)];
    x3 = [x3; x(3*(i-1) + 3)];
end


%plot(x1);
hold on;
grid on;
%plot(x2);
plot(x3);
plot(u);
title('u and x3 with r=15');