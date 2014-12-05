function lipmActuated(verify_manual_barrier_function)
% parameters
if nargin < 1
  verify_manual_barrier_function = false;
end
options.verify_manual_barrier_function = verify_manual_barrier_function;
options.plotfun = @visualize;

% path setup
addpath(fullfile('util'));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

% dynamics
nstates = 2;
f = @lipmDynamics;

% input limits
u_min = -1;
u_max = 1;

% failed states
x_ic_dist = 3;
g_Xfailed = @(x) (x(1) + x(2))^2 - x_ic_dist^2;

% initial states
x_star = [0; 0];
g_X0 = @(x) -(x - x_star)' * (x - x_star);

% manual barrier function
if verify_manual_barrier_function
  options.B_manual = @(x) (x(1) + x(2))^2 / (u_max)^2 - 1;
end

X0_margin = 1;

reset = [];
s_min = [];
s_max = [];
g_Xguard = [];
[B, u] = capturabilityBarrierRecursive(g_X0, f, nstates, u_min, u_max, reset, s_min, s_max, g_Xguard, g_Xfailed, X0_margin, options);

end



