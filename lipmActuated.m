function lipmActuated(N, verify_manual_barrier_function)
% parameters
if nargin < 1
  N = 1;
end
if nargin < 2
  verify_manual_barrier_function = true;
end
options.verify_manual_barrier_function = verify_manual_barrier_function;
options.plotfun = @visualize;

% path setup
addpath(fullfile('util'));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

% input limits
u_min = -1;
u_max = 1;

% failed states
x_ic_dist = 3;
g_Xfailed = @(x) (x(1) + x(2))^2 - x_ic_dist^2;

% initial states
x_star = [0; 0];
g_Xstar = @(x) -(x(1:2) - x_star)' * (x(1:2) - x_star);

% Discrete input limits
s_min = -1;
s_max = 1;

% Guard
t_min = 1;
g_Xguard = @(x) x(3) - t_min;

% manual barrier function

if verify_manual_barrier_function
%   dN = captureLimit(t_min, u_max, s_max, N);
  dN = captureLimit(t_min, u_max, s_max, N - 1);
  options.B_manual = @(x) (x(1) + x(2))^2 / (dN)^2 - 1;
  options.s_manual = @(x) -s_max * (x(1) + x(2)) / dN;
end

if N > 0
  dN_minus_one = captureLimit(t_min, u_max, s_max, N - 1);
  BN_minus_one_manual = @(x) (x(1) + x(2))^2 / (dN_minus_one)^2 - 1;
  f = @(x, u) [lipmDynamics(x, u); 1];
  nstates = 3;
  reset = @(x, s) [x(1) + s; x(2); 0];
else
  BN_minus_one_manual = [];
  f = @lipmDynamics;
  nstates = 2;
  reset = [];
end

[BN, uN] = capturabilityBarrier(BN_minus_one_manual, f, nstates, u_min, u_max, reset, s_min, s_max, g_Xguard, g_Xfailed, g_Xstar, options);

end

