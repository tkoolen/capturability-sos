function [B_fun, u_fun] = viableCaptureContinuous(f, nstates, u_min, u_max, g_Xtarget, g_Xfailed, g_Xstar, options)

% options
solver_options = spot_sdp_default_options();
solver_options.do_facial_reduction = false;
solver_options.verbose = 0;

bilinear_solve_options.max_iters = 150; %75;
bilinear_solve_options.rank_tol = 1e-5;

use_stored_initial_guess = true;
verify_manual_barrier_function = isfield(options, 'B_manual');
barrier_grow_iters = 15;

% spotless setup
solver = getSolver();
prog = spotsosprog;
bilinear_sos_constraints = cell(0);

% Indeterminates
[prog, x] = prog.newIndeterminate('x', nstates);

% Fill in indeterminates
g_Xtarget = functionHandleSubs(g_Xtarget, x);
g_Xfailed = functionHandleSubs(g_Xfailed, x);
g_Xstar = functionHandleSubs(g_Xstar, x);

% Inputs
u_degree = 1;
ninputs = length(u_min);
[prog, u] = prog.newFreePoly(monomials(x, 0 : u_degree), ninputs);

% Barrier function
B_degree = 2;
if verify_manual_barrier_function
  B = options.B_manual(x);
  [prog, B_full] = prog.newFreePoly(monomials(x, 0 : B_degree));
  prog = prog.withPolyEqs(B - B_full);
else
  [prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));
end

% Unsafe set constraint
[prog, lambda_Xfailed] = prog.newSOSPoly(monomials(x, 0 : 2));
prog = prog.withSOS(B - lambda_Xfailed * g_Xfailed - 1e-7); % B > 0 on Xf

% Barrier function derivative constraint
Bdot = diff(B, x) * f(x, u);
nu_B_degree = 2;
lambda_Xg_degree = 2;
[prog, nu_B] = prog.newFreePoly(monomials(x, 0 : nu_B_degree));
if ~isempty(g_Xtarget)
  [prog, z] = prog.newFreePoly(monomials(x, 0 : 4));
  for i = 1 : length(g_Xtarget)
    [prog, lambda_Xg] = prog.newSOSPoly(monomials(x, 0 : lambda_Xg_degree));
    prog = prog.withSOS(z - lambda_Xg * g_Xtarget(i)); % z >= 0 wherever g_Xtarget(i) >= 0
  end
  [prog, lambda_z] = prog.newSOSPoly(monomials(x, 0 : 2));
  bilinear_sos_constraints{end + 1} = -Bdot + nu_B * B + lambda_z * z + 1e-7; % Bdot(x) < 0 wherever B(x) = 0 and z(x) < 0
else
  bilinear_sos_constraints{end + 1} = -Bdot + nu_B * B + 1e-7; % Bdot(x) < 0 wherever B(x) = 0
end

% Input limits
[prog, nu_u_min_B] = prog.newFreePoly(monomials(x, 0 : nu_B_degree));
[prog, nu_u_max_B] = prog.newFreePoly(monomials(x, 0 : nu_B_degree));
if ~isempty(g_Xtarget)
  [prog, lambda_u_min_z] = prog.newSOSPoly(monomials(x, 0 : 2));
  bilinear_sos_constraints{end + 1} = u - u_min + nu_u_min_B * B + lambda_u_min_z * z; % u >= u_min wherever B(x) = 0 and z(x) < 0
  [prog, lambda_u_max_z] = prog.newSOSPoly(monomials(x, 0 : 2));
  bilinear_sos_constraints{end + 1} = u_max - u + nu_u_max_B * B + lambda_u_max_z * z; % u <= u_max wherever B(x) = 0 and z(x) < 0
else
  bilinear_sos_constraints{end + 1} = u - u_min + nu_u_min_B * B; % u >= u_min wherever B(x) = 0
  bilinear_sos_constraints{end + 1} = u_max - u + nu_u_max_B * B; % u <= u_max wherever B(x) = 0
end

% % Initial condition constraint
[prog, L0] = prog.newSOSPoly(monomials(x, 0 : 2));
prog = prog.withSOS(-B - L0 * g_Xstar - 1); % B <= -X0_margin on Xstar

indeterminates = x;

% Solve
B_fun = [];
u_fun = [];
if verify_manual_barrier_function
  for i = 1 : length(bilinear_sos_constraints)
    % SOS constraints that would normally be bilinear are now just linear
    % since B is given
    prog = prog.withSOS(bilinear_sos_constraints{i});
  end
  sol = prog.minimize(0, solver, solver_options);
  disp(char(sol.status));
 
  B_sol = sol.eval(B);
  u_sol = sol.eval(u);
  options.visualizer.visualize(B_sol, x, u_sol, f);

  % TODO: if success...
  B_fun = makeFunction(B_sol, x);
  u_fun = makeFunction(u_sol, x);
  
  vars = full(double(sol.eval(prog.variables)));
  save 'initial_guess.mat' vars;
else
  % load initial guess
  if use_stored_initial_guess
    load('initial_guess.mat');
    bilinear_solve_options.initial_guess = vars + 0.2 * randn(size(vars));
  end
  
  % iteratively grow zero level set of barrier function
  prog_base = prog;
  for i = 1 : barrier_grow_iters
    [sol, success] = solveBilinear(prog, bilinear_sos_constraints, indeterminates, solver, solver_options, bilinear_solve_options);
    disp(['status: ' char(sol.status)]);
    disp(['success: ' num2str(success)])
    fprintf('\n\n');
    
    B_sol = sol.eval(B);
    u_sol = sol.eval(u);
    
    if success
      vars_best = full(double(sol.eval(prog.variables)));
      options.visualizer.visualize(B_sol, x, u_sol, f);
      
      B_fun = makeFunction(B_sol, x);
      u_fun = makeFunction(u_sol, x);

      vars = full(double(sol.eval(prog.variables)));
      save 'initial_guess.mat' vars;

      prog = prog_base;
      [prog, L0] = prog.newSOSPoly(monomials(x, 0 : 2));
      prog = prog.withSOS(-B - L0 * -B_sol);
      
      bilinear_solve_options.initial_guess = vars_best + 0.1 * randn(size(vars_best));
    end
  end
end
end

