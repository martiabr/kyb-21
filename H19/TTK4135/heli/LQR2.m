%Weight and gain matrices for LQR
R_lqr = diag([20; 1]);
Q_lqr = diag([200; 20; 10; 5; 500; 10]);
[K_lqr,S,e] = dlqr(A_d,B_d,Q_lqr,R_lqr);    