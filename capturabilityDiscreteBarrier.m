function [B_fun, s_fun] = capturabilityDiscreteBarrier(B_prev, nstates, reset, s_min, s_max, g_Xguard, g_Xfailed, g_Xstar, options)

% options
solver_options = spot_sdp_default_options();
solver_options.do_facial_reduction = false;
solver_options.verbose = 0;

bilinear_solve_options.max_iters = 150; %75;
bilinear_solve_options.rank_tol = 1e-5;

use_stored_initial_guess = true;
verify_manual_barrier_function = isfield(options, 'B_manual');
barrier_grow_iters = 5;

% degrees
B_degree = 2;
s_degree = 1;

% spotless setup
solver = getSolver();
prog = spotsosprog;
bilinear_sos_constraints = cell(0);

% Indeterminates
[prog, x] = prog.newIndeterminate('x', nstates);

% Barrier function
if verify_manual_barrier_function
  B = options.B_manual(x);
  [prog, B_full] = prog.newFreePoly(monomials(x, 0 : B_degree));
  prog = prog.withPolyEqs(B - B_full);
else
  [prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));
end

% Unsafe set constraint
[prog, LF] = prog.newSOSPoly(monomials(x, 0 : 2));
prog = prog.withSOS(B - LF * g_Xfailed(x));

% discrete input
ns = length(s_min);
if verify_manual_barrier_function
  s = options.s_manual(x);
  [prog, s_full] = prog.newFreePoly(monomials(x, 0 : s_degree), ns);
  prog = prog.withPolyEqs(s - s_full);
else
  [prog, s] = prog.newFreePoly(monomials(x, 0 : s_degree), ns);
end

% discrete input limits
[prog, Ns_min] = prog.newFreePoly(monomials(x, 0 : 2));
[prog, Ls_guard_min] = prog.newSOSPoly(monomials(x, 0 : 2));
bilinear_sos_constraints{end + 1} = s - s_min + Ns_min * B - Ls_guard_min * g_Xguard(x); % s >= s_min on B(x) = 0, g_Xguard(x) >= 0
[prog, Ns_max] = prog.newFreePoly(monomials(x, 0 : 2));
[prog, Ls_guard_max] = prog.newSOSPoly(monomials(x, 0 : 2));
bilinear_sos_constraints{end + 1} = s_max - s + Ns_max * B - Ls_guard_max * g_Xguard(x); % s <= s_max on B(x) = 0, g_Xguard(x) >= 0

% B_prev(xPrime) <= 0 wherever xPrime = reset(x, s), B(x) <= 0, g_Xguard <= 0
[prog, x_prime] = prog.newIndeterminate('y', nstates);
indeterminates = [x; x_prime];

[prog, LB_prev_B] = prog.newFreePoly(monomials(indeterminates, 0 : 2));
[prog, LB_prev_guard] = prog.newSOSPoly(monomials(indeterminates, 0 : 2));
[prog, NB_prev_reset] = prog.newFreePoly(monomials(indeterminates, 0 : 2), nstates);
bilinear_sos_constraints{end + 1} = -B_prev(x_prime) + LB_prev_B * B - NB_prev_reset' * (x_prime - reset(x, s)) - LB_prev_guard * g_Xguard(x);

[prog, LB_prev] = prog.newSOSPoly(monomials(x, 0 : 2));
[prog, LB_prev_restrict] = prog.newSOSPoly(monomials(x, 0 : 2));
prog = prog.withSOS(-B - LB_prev * -B_prev(x) - LB_prev_restrict * g_Xguard(x)); % B 0-sublevel set is at least as large as B_prev 0-sublevel set restricted to guard set

% Solve
B_fun = [];
s_fun = [];
if verify_manual_barrier_function
  for i = 1 : length(bilinear_sos_constraints)
    % SOS constraints that would normally be bilinear are now just linear
    % since B is given
    prog = prog.withSOS(bilinear_sos_constraints{i});
  end
  sol = prog.minimize(0, solver, solver_options);
  disp(char(sol.status));
 
  B_sol = sol.eval(B);
  s_sol = sol.eval(s);
  options.plotfun(B_sol, x, s_sol, reset);

  % TODO: if success...
  B_fun = makeFunction(B_sol, x);
  s_fun = makeFunction(s_sol, x);
  
  vars = full(double(sol.eval(prog.variables)));
  save 'initial_guess.mat' vars;
else
  % load initial guess
  if use_stored_initial_guess
    load('initial_guess.mat');
    bilinear_solve_options.initial_guess = vars; % + 0.5 * randn(size(vars));
  end
  
  % iteratively grow zero level set of barrier function
  prog_base = prog;
  for i = 1 : barrier_grow_iters
    [sol, success] = solveBilinear(prog, bilinear_sos_constraints, indeterminates, solver, solver_options, bilinear_solve_options);
    disp(['status: ' char(sol.status)]);
    disp(['success: ' num2str(success)])
    fprintf('\n\n');
    
    B_sol = clean(sol.eval(B), 1e-5);
    s_sol = clean(sol.eval(s), 1e-5);
    
    if success
      vars_best = full(double(sol.eval(prog.variables)));
      options.plotfun(B_sol, x, s_sol, reset);
      
      B_fun = makeFunction(B_sol, x);
      s_fun = makeFunction(s_sol, x);
      
      prog = prog_base;
      [prog, L0] = prog.newSOSPoly(monomials(x, 0 : 2));
      prog = prog.withSOS(-B - L0 * -B_sol);
      
      bilinear_solve_options.initial_guess = vars_best + 0.1 * randn(size(vars_best));
    end
  end
end

end

