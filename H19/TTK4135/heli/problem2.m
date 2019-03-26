%% Initialization and model definition
init03;

% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A1 = [1 delta_t 0 0; 
    0 1 -delta_t*K_2 0; 
    0 0 1 delta_t; 
    0 0 -delta_t*K_1*K_pp 1-(delta_t*K_1*K_pd)];

B1 = [0; 0; 0; delta_t*K_1*K_pp];

% Number of states and inputs
mx = size(A1,2);
mu = size(B1,2);

% Initial values
x1_0 = pi;
x2_0 = 0;
x3_0 = 0;
x4_0 = 0;
x0 = [x1_0 x2_0 x3_0 x4_0]';

% Time horizon and initialization
N  = 100;
M  = N;
z  = zeros(N*mx+M*mu,1);
z0 = z;

% Bounds
ul 	    = -30*pi/180;
uu 	    = 30*pi/180;

xl      = -Inf*ones(mx,1);
xu      = Inf*ones(mx,1);
xl(3)   = ul;
xu(3)   = uu;

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu);
vlb(N*mx+M*mu)  = 0;
vub(N*mx+M*mu)  = 0;

% Generate the matrix Q and the vector c
Q1 = zeros(mx,mx);
Q1(1,1) = 1;
Q1(2,2) = 0;
Q1(3,3) = 0;
Q1(4,4) = 0;
P1 = 10;
Q = gen_q(Q1,P1,N,M);
c = zeros(500, 1);

%% Generate system matrixes for linear model
Aeq = gen_aeq(A1,B1,N,mx,mu);
beq = [A1*x0; zeros(396,1)];

%% Solve QP problem with linear model
tic
[z,lambda] = quadprog(Q, c, [], [], Aeq, beq, ...
vlb, vub);
t1 = toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)];

x1 = [x0(1);z(1:mx:N*mx)];
x2 = [x0(2);z(2:mx:N*mx)];
x3 = [x0(3);z(3:mx:N*mx)];
x4 = [x0(4);z(4:mx:N*mx)];

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
