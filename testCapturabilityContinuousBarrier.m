function testCapturabilityContinuousBarrier(N, verify_manual_barrier_function)

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
s_max = 1;

% dynamics
f = @(x, u) [lipmDynamics(x, u); 1];

% input limits
u_min = -1;
u_max = 1;

if N > 0
  dN_minus_one = captureLimit(t_min, u_max, s_max, N - 1);
  B_prev = @(x) (x(1) + x(2))^2 / (dN_minus_one + s_max)^2 - 1;
else
  B_prev = @(x) 0;
end

nstates = 3;

% failed states
x_ic_dist = 3;
g_Xf = @(x) (x(1) + x(2))^2 - x_ic_dist^2;

options.plotfun = @visualizeLIPM;

if verify_manual_barrier_function
  if N == 0
    dN = captureLimit(t_min, u_max, s_max, N);
    options.B_manual = @(x) (x(1) + x(2))^2 / dN^2 - 1;
  else
    %   options.B_manual = B_prev;
    
    %   dN_minus_one = captureLimit(t_min, u_max, s_max, N - 1);
    %   options.B_manual = @(x) (x(1) + x(2))^2 / (dN_minus_one + s_max)^2 - 1;
    %   options.s_manual = @(x) -s_max * (x(1) + x(2)) / (dN_minus_one + s_max);
  end
end

g_Xg = @(x) -B_prev(x);
[B_fun, u_fun] = viableCaptureContinuous(f, nstates, u_min, u_max, g_Xg, g_Xf, options);

end