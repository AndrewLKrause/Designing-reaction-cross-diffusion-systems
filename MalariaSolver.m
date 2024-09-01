function [U,x,ui,vi] = ...
    MalariaSolver(dims, m, params, tols, T,showProgBar)
% This code solves the 3-component cross-diffusion malaria system on a
% square domain in 1D and 2D.

% Parameters of the model.
[L, b, b_H, d_H, b_M, d_M, c, r, Q, D, eta, BC] = deal(params{:});
%L,  b, b_H, d_H, b_M, d_M, c, r, Q, D, eta,  BC

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
%    'warning: this model can behave badly in 2D!'
end

% Indices corresponding to u variable and v variable. THis lets us stack
% them both in a vector U and write u = U(ui) and v = U(vi).
ui = 1:N; vi = N+1:2*N; wi = 2*N+1:3*N;

% Reaction kinetics & nonlinear diffusion coefficient

F = @(u,v,w) (b_H-d_H)*u-c*w.*u+r*v;
G = @(u,v,w) -d_H*v + c*w.*u -r*v;
H = @(u,v,w) -d_M*w+b*(Q-w).*v;
xi = @(u)2*u./(1+u.^2);

% Homogeneous initial conditions, needed to scale the nonlinear diffusion.

Hss  =  d_H*d_M*(d_H+r)/(b*(d_H*(c*Q+d_H+r)-b_H*(d_H+r)));
Iss  =  Hss*(b_H-d_H)/d_H;
Pss  =  (b_H-d_H)*(d_H+r)/(c*d_H);



if (dims == 1)
    % 1D Laplacian
    Adv = (1/(2*dx))*Adv;
    Lap = (1/dx)^2*Lap;

    NLD = @(Du,U)Du.*(Lap*U) + (Adv*Du.*(Adv*U));

elseif (dims == 2)
    % 2D Laplacian
    I = speye(m);
    Advx = (1/(2*dx))*kron(Adv,I);
    Advy = (1/(2*dx))*kron(I, Adv);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));

    NLD = @(Du,U)Du.*(Lap*U) + (Advx*Du).*(Advx*U) + (Advy*Du).*(Advy*U);

end


F = @(t,U)[F(U(ui),U(vi),U(wi)) + D(1,1)*Lap*U(ui)+D(2,1)*NLD(xi(U(ui)./Hss),U(vi))+D(3,1)*NLD(xi(U(ui)./Hss),U(wi));
    G(U(ui),U(vi),U(wi)) + D(1,2)*NLD(xi(U(vi)./Iss),U(ui))+D(2,2)*Lap*U(vi)+D(3,2)*Lap*U(wi);
    H(U(ui),U(vi),U(wi)) + D(1,3)*Lap*U(ui)+D(2,3)*Lap*U(vi)+D(3,3)*Lap*U(wi)];


% Initial condition - this is a small normally distributed perturbation of
% the homogeneous steady state of our kinetics
U0 = [Hss*(1 + eta*randn(N,1)); Iss*(1 + eta*randn(N,1)); Pss*(1 + eta*randn(N,1))];

% This is the Jacobian sparsity pattern. That is, if you compute the
% Jacobian of the vector function F above for the vector argument U, this
% matrix is where all of the nonzero elements are. This is important for
% implicit timestepping!
I = speye(N); Z = sparse(zeros(N));
%Warning: I have used D' in writing out the elements above, hence why this
%pattern uses that for the coupling.
JacSparse = sparse([I, I, I; I, I, I; Z, I, I]+abs(kron(D',Lap)));
odeOptions = odeset('JPattern',JacSparse,'RelTol',tols,'AbsTol',tols,'InitialStep',1e-6);
if (showProgBar)
    odeOptions = odeset(odeOptions,'OutputFcn',@ODEProgBar);
end

% Solve the system using an implicit stiff timestepper.
[~,U] = ode15s(F,T,U0,odeOptions);
end


%Nonlinear