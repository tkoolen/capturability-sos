function testTrivialSystem(N, verify_manual_barrier_function)
% parameters
if nargin < 1
  N = 0;
end
if nargin < 2
  verify_manual_barrier_function = false;
end
options.verify_manual_barrier_function = verify_manual_barrier_function;
options.plotfun = @visualizeTrivialSystem;

% path setup
addpath(fullfile('util'));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

% input limits
u_min = 0;
u_max = 0;

% failed states
g_Xfailed = @(x) x^2 - 3^2;

% initial states
g_Xstar = @(x) -x^2;

% manual barrier function
B0_manual = @(x) x^2 - 1^2;


% Discrete input limits
s_min = -1;
s_max = 1;

% Guard
g_Xguard = @(x) 1;

nstates = 1;

% zero-step capturability
if N == 0
  % dynamics
  f = @(x, u) 0;
  reset = [];
  if verify_manual_barrier_function
    options.B_manual = B0_manual;
  end
  [B, u] = capturabilityBarrier([], f, nstates, u_min, u_max, [], s_min, s_max, g_Xguard, g_Xfailed, g_Xstar, options);
else % 1-step capturability test
  % dynamics
  f = @(x, u) 0;
  
  % Reset map
  reset = @(x, s) x + s;

  if verify_manual_barrier_function
    B1_manual = @(x) x^2 - 1.5^2;
    options.B_manual = B1_manual;
    options.s_manual = @(x) -x / 2;
  end
  
  [B, u, s] = capturabilityBarrier(B0_manual, f, nstates, u_min, u_max, reset, s_min, s_max, g_Xguard, g_Xfailed, g_Xstar, options);
end

end
