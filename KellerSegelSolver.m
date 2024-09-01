function [U,x,ui,vi] = ...
    KellerSegelSolver(dims, m, params, tols, T,showProgBar)
% This code solves the Keller-Segel chemotaxis system on a square
% domain in 1D and 2D.

% Parameters of the model.
[L, a, b, c, d, e, f, g, h,i, epsilon, eta, BC] = deal(params{:});

% Spatial step size
dx = L/(m-1);

% Spatial domain (needed for plotting only)
x = linspace(0,L,m);

% Useful sparse matrices for transport
eVec = ones(m,1);
Lap = spdiags([eVec,-2*eVec,eVec],[1,0,-1],m,m);
Adv = spdiags([eVec,-eVec],[1,-1],m,m);

if(BC==0)
    % Neumann boundary conditions
    Lap(1,1) = -1; Lap(end,end) = -1;
    Adv(1,2)=0; Adv(end,end-1)=0;
else
    % Periodic boundary conditions
    Lap(1,end) = 1; Lap(end,1) = 1;
    Adv(end,1) = 1; Adv(1,end) = -1;
end
if(dims==1)
    N = m;
elseif(dims==2)
    N=m^2;
    'warning: this model can behave badly in 2D!'
end

% Indices corresponding to u variable and v variable. THis lets us stack
% them both in a vector U and write u = U(ui) and v = U(vi).
ui = 1:N; vi = N+1:2*N; wi = 2*N+1:3*N;

% Reaction kinetics & nonlinear diffusion coefficient
F = @(u,v,w) b.*u.*(1 - u);
G = @(u,v,w) c*u + e*v+f*w-epsilon*v.^3;
H = @(u,v,w) g*u + h*v+i*w-epsilon*w.^3;
xi = @(u)2*a*u./(1+u.^2);

% Set the Laplacian and the nonlinear diffusion DivGrad based on dimension.
% Write the full RHS based on dimension.
if (dims == 1)
    % 1D Laplacian
    Adv = (1/(2*dx))*Adv;
    Lap = (1/dx)^2*Lap;

    NLD = @(D,U)D.*(Lap*U) + (Adv*D.*(Adv*U));

elseif (dims == 2)
    % 2D Laplacian
    I = speye(m);
    Advx = (1/(2*dx))*kron(Adv,I);
    Advy = (1/(2*dx))*kron(I, Adv);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));

    NLD = @(D,U)D.*(U) + (Advx*D).*(Advx*U) + (Advy*D).*(Advy*U);

end

F = @(t,U)[F(U(ui),U(vi),U(wi)) + d*Lap*U(ui) - a*NLD(xi(U(ui)),U(vi));
    G(U(ui),U(vi),U(wi)) + Lap*U(vi);
    H(U(ui),U(vi),U(wi)) + Lap*U(wi)];

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
U0 = [1 + eta*randn(N,1); 0 + eta*randn(N,1); 0 + eta*randn(N,1)];

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
I = speye(N); Z = sparse(zeros(N));
JacSparse = sparse([Lap, Lap, Z; I, Lap, I; I, I, Lap]);
odeOptions = odeset('JPattern',JacSparse,'RelTol',tols,'AbsTol',tols,'InitialStep',1e-6);
if (showProgBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[~,U] = ode15s(F,T,U0,odeOptions);
end