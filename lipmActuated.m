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
g_Xstar = @(x) -(x - x_star)' * (x - x_star);

% manual barrier function
B_manual = @(x) (x(1) + x(2))^2 / (u_max)^2 - 1;
if verify_manual_barrier_function
  options.B_manual = B_manual;
end

% Reset map
reset = @(x, s) [x(1) + s; x(2); 0];

% Discrete input limits
s_min = -1;
s_max = 1;

% Guard
t_min = 1;
g_Xguard = @(x) -x(3) + t_min;

% zero-step capturability
zero_step = true;
if zero_step
  [B, u] = capturabilityBarrier([], f, nstates, u_min, u_max, [], s_min, s_max, g_Xguard, g_Xfailed, g_Xstar, options);
else
  % 1-step capturability test
  [B, u] = capturabilityBarrier(B_manual, f, nstates, u_min, u_max, reset, s_min, s_max, g_Xguard, g_Xfailed, g_Xstar, options);
end


end



