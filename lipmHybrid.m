function lipmHybrid()

% path setup
addpath(fullfile('util'));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

options.solver_options = spot_sdp_default_options();
options.solver_options.do_facial_reduction = false;
options.solver_options.verbose = 0;
options.bilinear_solve_options.max_iters = 75;
options.bilinear_solve_options.rank_tol = 1e-3;
options.bilinear_solve_options.psd_constraint_size = 2;
options.solver = getSolver();

% Dynamics
nstates = 3;
f = @(x, u) [lipmDynamics(x(1 : 2), u); 1];

% continuous input limits
u_min = -1;
u_max = 1;

% Reset map
reset = @(x, s) [x(1) + s; x(2); 0];

% Discrete input limits
s_min = -1;
s_max = 1;

% Guard
t_min = 1;
g_Xguard = @(x) -x(3) + t_min;

% Failed states
x_ic_dist = 3;
g_Xfailed = @(x) (x(1) + x(2))^2 - x_ic_dist^2;

% Barrier function for previous level
B_prev = @(x) (x(1) + x(2))^2 / (u_max)^2 - 1;

B = capturabilityBarrier(B_prev, f, nstates, u_min, u_max, reset, s_min, s_max, g_Xguard, g_Xfailed, options);

end