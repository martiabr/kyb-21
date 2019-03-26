%% Initialization and model definition
init03;

% Discrete time system model. x = [lambda r p p_dot e e_dot]'
T	= 0.25; % sampling time
A_c = [0 1 0 0 0 0; 
      0 0 -K_2 0 0 0; 
      0 0 0 1 0 0; 
      0 0 -K_1*K_pp -K_1*K_pd 0 0;
      0 0 0 0 0 1;
      0 0 0 0 -K_3*K_ep -K_3*K_ed];

B_c = [0 0; 0 0; 0 0; K_1*K_pp 0; 0 0; 0 K_3*K_ep];

A_d = eye(6) + T*A_c;
B_d = T*B_c;

% Number of states and inputs
mx = size(A_d,2);
mu = size(B_d,2);

% Initial values
x1_0 = pi;
x2_0 = 0;
x3_0 = 0;
x4_0 = 0;
x5_0 = 0;
x6_0 = 0;
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';

% Time horizon and initialization
N  = 60;
M  = N;
z  = zeros(N*mx+M*mu,1);
z0 = z;

% Bounds
ul      = -Inf*ones(mu,1);
uu      = Inf*ones(mu,1); 
ul(1)   = -0.45;
uu(1) 	= 0.45;
ul(2)   = -0.3;
uu(2)   = 0.3;

xl      = -Inf*ones(mx,1);
xu      = Inf*ones(mx,1);
xl(3)   = -0.45;
xu(3)   = 0.45;
xl(5)   = -0.3;
xu(5)   = 0.3;
xl(6)   = -0.2;
xu(6)   = 0.2;

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu);
vlb(N*mx+M*mu)  = 0;
vub(N*mx+M*mu)  = 0;

% Generate the matrix Q and the vector c
Q1 = zeros(mx);
Q1(1,1) = 1;
Q1(2,2) = 0;
Q1(3,3) = 0;
Q1(4,4) = 0;
Q1(5,5) = 0;
Q2(6,6) = 0;

P1 = zeros(mu);                     
P1(1,1) = 0.1;
P1(2,2) = 0.1;

Q = gen_q(Q1,P1,N,M);
c = zeros(N*(mu+mx), 1);

%% Generate system matrixes for linear model
Aeq = gen_aeq(A_d,B_d,N,mx,mu);
beq = [A_d*x0; zeros((N-1)*mx,1)];

%% Solve nonlinear problem
opts = optimset('Display','iter','Algorithm','sqp',...
'MaxFunEval',inf,'MaxIter',Inf);

f = @(z) 0.5*z'*Q*z;

[z, lambda] = fmincon(f, z0, [], [], Aeq, beq, ...
vlb, vub, @nonlcon, opts); 

%% Extract control inputs and states
u1 = z(N*mx+1:mu:N*mx+M*mu);
u2 = z(N*mx+2:mu:N*mx+M*mu);

x1 = [x0(1);z(1:mx:(N-1)*mx)];
x2 = [x0(2);z(2:mx:(N-1)*mx)];
x3 = [x0(3);z(3:mx:(N-1)*mx)];
x4 = [x0(4);z(4:mx:(N-1)*mx)];
x5 = [x0(5);z(5:mx:(N-1)*mx)];
x6 = [x0(6);z(6:mx:(N-1)*mx)];

padding = 8/T;  %used 5 here
zero_padding = zeros(padding,1);
unit_padding  = ones(padding,1);

u1   = [zero_padding; u1; zero_padding];
u2   = [zero_padding; u2; zero_padding];

x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];