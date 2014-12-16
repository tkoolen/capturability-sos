function testCapturabilityDiscreteBarrier(N, verify_manual_barrier_function)

% path setup
addpath(fullfile('util'));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

% parameters
if nargin < 1
  N = 1;
end
if nargin < 2
  verify_manual_barrier_function = true;
end

t_min = 1;

% Discrete input limits
s_max = 1;
s_min = -s_max;

u_max = 1;

dN = captureLimit(t_min, u_max, s_max, N);
B_prev = @(x) (x(1) + x(2))^2 / (dN)^2 - 1;
nstates = 3;

% reset map
reset = @(x, s) [x(1) + s; x(2); 0];

% Guard
g_Xguard = @(x) x(3) - t_min;

% failed states
x_ic_dist = 3;
g_Xfailed = @(x) (x(1) + x(2))^2 - x_ic_dist^2;

g_Xstar = [];

options = struct;

options.plotfun = @(varargin) [];

if verify_manual_barrier_function
%   options.B_manual = @(x) B_prev(x);
%   options.s_manual = @(x) 0;
  
  dN_minus_one = captureLimit(t_min, u_max, s_max, N - 1);
  options.B_manual = @(x) (x(1) + x(2))^2 / (dN_minus_one + s_max)^2 - 1;
  options.s_manual = @(x) -s_max * (x(1) + x(2)) / (dN_minus_one + s_max);
end

[B_fun, s_fun] = capturabilityDiscreteBarrier(B_prev, nstates, reset, s_min, s_max, g_Xguard, g_Xfailed, g_Xstar, options);

end