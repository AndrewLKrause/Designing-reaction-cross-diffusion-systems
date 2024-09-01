function [U,x,ui,vi] = ...
    HyperbolicSolver(dims, m, params, tols, T,showProgBar)
% This code solves the 3-component cross-diffusion Schnakenberg system on a
% square domain in 1D and 2D.

% Parameters of the model.
[L, a, b, delta, d11, d12, d21, d22, tau, eta, BC] = deal(params{:});



% Spatial step size
dx = L/(m-1);

% Spatial domain (needed for plotting only)
x = linspace(0,L,m);

% (Sparse) Laplacian matrix
eVec = ones(m,1);
Lap = spdiags([eVec,-2*eVec,eVec],[1,0,-1],m,m);

if(BC==0)
    % Neumann boundary conditions
    Lap(1,1) = -1; Lap(end,end) = -1;
else
    % Periodic boundary conditions
    Lap(1,end) = 1; Lap(end,1) = 1;
end
if(dims==1)
    N = m;
elseif(dims==2)
    N=m^2;
end
% Indices corresponding to u variable and v variable. THis lets us stack
% them both in a vector U and write u = U(ui) and v = U(vi).
ui = 1:N; vi = N+1:2*N; wi = 2*N+1:3*N; qi = 3*N+1:4*N;

% Reaction kinetics
F = @(u,v,w,q) w;
G = @(u,v,w,q) q;
H = @(u,v,w,q) (u-a*u.^3+v-b-w)/tau;
J = @(u,v,w,q) (b-u-3*a*v-q)/tau;

% Set the Laplacian based on dimension.
if (dims == 1)
    % 1D Laplacian
    Lap = (1/dx)^2*Lap;
elseif (dims == 2)
    % 2D Laplacian
    I = speye(m);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));

end

%Write the RHS.
RHS = @(t,U)[F(U(ui),U(vi),U(wi),U(qi)) + delta*Lap*U(ui);
    G(U(ui),U(vi),U(wi),U(qi)) + delta*Lap*U(vi);
    H(U(ui),U(vi),U(wi),U(qi)) + (d11/tau)*Lap*U(ui)+(d12/tau)*Lap*U(vi)+Lap*U(wi);
    J(U(ui),U(vi),U(wi),U(qi)) + (d21/tau)*Lap*U(ui)+(d22/tau)*Lap*U(vi)+Lap*U(qi)];

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
Nss  =  (1  -  3*a)^2*a^6*(3*a*(27*a*b^2  -  4)+4);
Mss  =  (9*(1-3*a)*a^4*b+sqrt(Nss))^(1/3);
uss  =  ((2^(1/3))*Mss^2  +  6*a^3  -  2*a^2)/(3*2^(2/3)*a^2*Mss);
vss  =  (b-uss)/(3*a);

U0 = [uss*(1 + eta*randn(N,1)); vss*(1 + eta*randn(N,1)); 0*(1 + eta*randn(N,1));0*(1 + eta*randn(N,1))];



% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
I = speye(N); Z = sparse(zeros(N));
JacSparse = sparse([Lap, Z, I, Z; Z, Lap, Z, I; Lap, Lap, Lap, Z; Lap, Lap, Z, Lap]);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tols,'AbsTol',tols,'InitialStep',1e-6);
if (showProgBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[~,U] = ode15s(RHS,T,U0,odeOptions);
end