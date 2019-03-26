function [c, ceq] = nonlcon(z)

alpha = 0.2; 
beta = 20; 
lamda_t = 2*pi/3; 
N = 40;

c = zeros(N, 1);
ceq = [];

for i=1:N
    c(i) = alpha*exp(-beta*(z((i-1)*6 + 1) - lamda_t)^2)...
    - z((i-1)*6 + 5);
end