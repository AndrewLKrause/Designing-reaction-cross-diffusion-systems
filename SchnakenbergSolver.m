function [U,x,ui,vi] = ...
    SchnakenbergSolver(dims, m, params, tols, T,showProgBar)
% This code solves the 3-component cross-diffusion Schnakenberg system on a
% square domain in 1D and 2D.

% Parameters of the model.
[L, a, b, D, eta, BC] = deal(params{:});

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
ui = 1:N; vi = N+1:2*N; wi = 2*N+1:3*N;

% Reaction kinetics
F = @(u,v,w) a-u+u.^2.*v;
G = @(u,v,w) b-u.^2.*v;
H = @(u,v,w) u-w;

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
F = @(t,U)[F(U(ui),U(vi),U(wi)) + D(1,1)*Lap*U(ui)+D(2,1)*Lap*U(vi)+D(3,1)*Lap*U(wi);
    G(U(ui),U(vi),U(wi)) + D(1,2)*Lap*U(ui)+D(2,2)*Lap*U(vi)+D(3,2)*Lap*U(wi);
    H(U(ui),U(vi),U(wi)) + D(1,3)*Lap*U(ui)+D(2,3)*Lap*U(vi)+D(3,3)*Lap*U(wi)];

% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
U0 = [(a+b)*(1 + eta*randn(N,1)); (b/(a+b)^2)*(1 + eta*randn(N,1)); (a+b)*(1 + eta*randn(N,1))];

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
I = speye(N); Z = sparse(zeros(N));
JacSparse = sparse([I, I, Z; I, I, Z; I, Z, I]+abs(kron(D',Lap)));
odeOptions = odeset('JPattern',JacSparse,'RelTol',tols,'AbsTol',tols,'InitialStep',1e-6);
if (showProgBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[~,U] = ode15s(F,T,U0,odeOptions);
end