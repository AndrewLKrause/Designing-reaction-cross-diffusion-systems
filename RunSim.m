function [U,x,T,ui,vi] = RunSim(modelName,dims)

% Show a progress bar?
showProgBar = true;
% Set default random seed/
rng('default');

%Default variance of random seed
eta = 0.01;

% Number of gridpoints per dimension. Use 80-300 or so for 2D, and ideally
% 300-3000 or so for 1D depending on the structures that need to be
% resolved.
if(dims==1)
    m = 1000;
    number_dt_steps = 15000;
elseif(dims==2)
    m=200;
    number_dt_steps = 1000;
end

% Numerical tolerances (absolute and relative).
tols = 1e-9;

switch modelName
    case 'KellerSegelWaveNeumann'
        Params = {5, 3, 8, 33.75, 1.1, -512.5, -500, -0.41922, 512.577, 500, 0.01,    eta, 0};
        %        [L, a, b, c,     d,   e,      f,    g,        h,       i,   epsilon, eta,  BC (1 periodic, 0 Neumann)]
        Solver = @KellerSegelSolver;
        Tend=3;
    case 'KellerSegelWavePeriodic'
        Params = {5, 3, 8, 33.75, 1.1, -512.5, -500, -0.41922, 512.577, 500, 0.01,    eta, 1};
        %        [L, a, b, c,     d,   e,      f,    g,        h,       i,   epsilon, eta,  BC (1 periodic, 0 Neumann)]
        Solver = @KellerSegelSolver;
        Tend=3;
    case 'KellerSegelTuringPeriodic'
        Params = {10, 3, 1, 36.7, 0.9,   -32, -1, -1.5, -1618, -64, 0.01,    eta, 1};
        %        [L, a, b,  c,    d,     e,   f,  g,    h,     i,   epsilon, eta,  BC (1 periodic, 0 Neumann)]
        Solver = @KellerSegelSolver;
        Tend=3;
    case 'KellerSegelTuringNeumann'
        Params = {10, 3, 1, 36.7, 0.9,   -32, -1, -1.5, -1618, -64, 0.01,    eta, 0};
        %        [L, a, b,  c,    d,     e,   f,  g,    h,     i,   epsilon, eta,  BC (1 periodic, 0 Neumann)]
        Solver = @KellerSegelSolver;
        Tend = 15;
    case 'SchnakenbergWavePeriodic'
        delta = 0.8;
        Params = {2.38*pi, 1, 0.5, [delta, 0, 0; delta - 1, 1, 198 - 198*delta; 0, 0, delta]', eta, 1};
        %        [L,  a, b,   D,                              eta,  BC]
        Solver = @SchnakenbergSolver;
        Tend=60;
    case 'SchnakenbergWaveNeumann'
        delta = 0.8;
        Params = {2.38*pi, 1, 0.5, [delta, 0, 0; delta - 1, 1, 198 - 198*delta; 0, 0, delta]', eta, 0};
        %        [L,  a, b,   D,                              eta,  BC]
        Solver = @SchnakenbergSolver;
        Tend = 1e2;
    case 'MalariaTuringNeumann'
        Params = {6, 0.1, 1,   0.1, 1,   0.3, 0.25, 0.5, 100, [1, 0.5, 0.307225; 0.870348, 1, 0; 0, 0, 0.045]', eta, 0};
        %        [L,  b,   b_H, d_H, b_M, d_M, c,    r,    Q,   D,                                        eta,  BC]
        Solver = @MalariaSolver;
        Tend=100;
    case 'MalariaTuringPeriodic'
        Params = {6, 0.1, 1,   0.1, 1,   0.3, 0.25, 0.5, 100, [1, 0.5, 0.307225; 0.870348, 1, 0; 0, 0, 0.045]', eta, 1};
        %        [L,  b,   b_H, d_H, b_M, d_M, c,    r,    Q,   D,                                        eta,  BC]
        Solver = @MalariaSolver;
        Tend=100;
    case 'HyperbolicWavePeriodic'
        Params = {30, 0.257, 0.98, 1.3,    - 1,  - 1,  0,   - 2,  0.1,  eta, 1};
        %        [L,  a,     b,    delta,  d11, d12, d21, d22, tau,  eta,  BC]
        Solver = @HyperbolicSolver;
        Tend=150;
    case 'HyperbolicWaveNeumann'
        Params = {6, 0.257, 0.98, 1.3,    - 1,  - 1,  0,   - 2,  0.1,  eta, 0};
        %        [L,  a,     b,    delta,  d11, d12, d21, d22, tau,  eta,  BC]
        Solver = @HyperbolicSolver;
        Tend=150;

    otherwise
        disp('Unknown model.')
        return;
end


T = linspace(0,Tend,number_dt_steps);
[U,x,ui,vi] = Solver(dims, m, Params, tols, T,showProgBar);

%If everything is very flat, use the following to remove some of the
%initial perturbation from plotting to see smaller colour ranges.
%T = T(2:end);
%U = U(2:end,:);

end