function lipmActuated(verify_manual_barrier_function)
% path setup
addpath(fullfile('util'));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

% parameters
if nargin < 1
  verify_manual_barrier_function = false;
end
use_stored_initial_guess = true;
solver_options = spot_sdp_default_options();
solver_options.do_facial_reduction = false;
solver_options.verbose = 0;
bilinear_solve_options.max_iters = 75;
bilinear_solve_options.rank_tol = 1e-3;
u_degree = 2;
Nu_degree = 2;
B_degree = 2;
L0_degree = 2;
LF_degree = 2;
NBdot_degree = 2;

% spotless setup
solver = getSolver();
prog = spotsosprog;

% Indeterminates
nstates = 2;
ninputs = 1;
[prog, x] = prog.newIndeterminate('x', nstates);
r = x(1);
rd = x(2);

% Inputs
u_min = -1;
u_max = 1;
[prog, u] = prog.newFreePoly(monomials(x, 0 : u_degree), ninputs);

% Dynamics
f = @(x, u) lipmDynamics(x, u);

% Barrier function
if verify_manual_barrier_function
  B = (r + rd)^2 / (u_max)^2 - 1;
  [prog, B_full] = prog.newFreePoly(monomials(x, 0 : B_degree));
  prog = prog.withPolyEqs(B - B_full);
else
  [prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));
end

% Unsafe set constraint
[prog, LF] = prog.newSOSPoly(monomials(x, 0 : LF_degree));
% rf_dist = 10000;
% g_Xu = r' * r - rf_dist^2;
x_ic_dist = 3;
g_Xu = (r + rd)' * (r + rd) - x_ic_dist^2;
prog = prog.withSOS(B - LF * g_Xu); % B >= 0 on g_Xu

% Input limits
bilinear_sos_constraints = cell(0);
[prog, Nu_min] = prog.newFreePoly(monomials(x, 0 : Nu_degree));
bilinear_sos_constraints{end + 1} = u - u_min + Nu_min * B; % u >= u_min on B(x) = 0
[prog, Nu_max] = prog.newFreePoly(monomials(x, 0 : Nu_degree));
bilinear_sos_constraints{end + 1} = u_max - u + Nu_max * B; % u <= u_max on B(x) = 0

% Barrier function derivative constraint
Bdot = diff(B, x) * f(x, u);
[prog, NBdot] = prog.newFreePoly(monomials(x, 0 : NBdot_degree));
epsilon = 1e-10;
bilinear_sos_constraints{end + 1} = -Bdot + NBdot * B - epsilon; % Bdot <= -epsilon on B(x) = 0

% Initial condition constraint
prog_base = prog;
x_star = [0; 0];
g_X0 = - (x - x_star)' * (x - x_star);
% g_X0 = -(r + rd)' * (r + rd);
X0_margin = 1;
[prog, L0] = prog.newSOSPoly(monomials(x, 0 : L0_degree));
prog = prog.withSOS(-B - L0 * g_X0 - X0_margin); % B <= -X0_margin on g_X0

% Solve
if verify_manual_barrier_function
  for i = 1 : length(bilinear_sos_constraints)
    % SOS constraints that would normally be bilinear are now just linear
    % since B is given
    prog = prog.withSOS(bilinear_sos_constraints{i});
  end
  sol = prog.minimize(0, solver, solver_options);
  disp(char(sol.status));
  
  B_sol = sol.eval(B);
  Bdot_sol = sol.eval(Bdot);
  u_sol = sol.eval(u);
  visualize(B_sol, Bdot_sol, x, u_sol, f);
  
  % store initial guess
  initial_guess = {double(sol.eval([decomp(u, x); decomp(B_full, x); decomp(Nu_min, x); decomp(Nu_max, x); decomp(NBdot, x)]))};
  save 'initial_guess.mat' initial_guess;
else
  % load initial guess
  if use_stored_initial_guess
    load('initial_guess.mat');
    bilinear_solve_options.initial_guess = addNoise(initial_guess, 0.5);
  end
  
  % iteratively grow zero level set of barrier function
  for i = 1 : 20
    [sol, success, sol_w] = solveBilinear(prog, bilinear_sos_constraints, x, solver, solver_options, bilinear_solve_options);
    disp(['status: ' char(sol.status)]);
    disp(['success: ' num2str(success)])
    fprintf('\n\n');
    
    B_sol = sol.eval(B);
    Bdot_sol = sol.eval(Bdot);
    u_sol = sol.eval(u);
    
    if success
      sol_w_best = sol_w;
      X0_margin = 0;
      g_X0 = -B_sol;
      visualize(B_sol, Bdot_sol, x, u_sol, f);
    end
    
    prog = prog_base;
    [prog, L0] = prog.newSOSPoly(monomials(x, 0 : L0_degree));
    prog = prog.withSOS(-B - L0 * g_X0 - X0_margin);
    
    if exist('sol_w_best', 'var')
      bilinear_solve_options.initial_guess = addNoise(sol_w_best, 1);
    end
  end
end

end

function w_noisy = addNoise(w, sigma)
w_noisy = w;
for row = 1 : size(w, 1)
  for col = 1 : size(w, 2)
    w_noisy{row, col} = w_noisy{row, col} + sigma * randn;
  end
end
end


