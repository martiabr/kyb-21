%Weight and gain matrices for LQR
R_lqr = 1;
Q_lqr = diag([100; 1; 10; 1;]);
[K_lqr,S,e] = dlqr(A1,B1,Q_lqr,R_lqr);
