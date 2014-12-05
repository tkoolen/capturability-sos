function [B_fun, u_fun] = capturabilityBarrier(B_prev, f, nstates, u_min, u_max, reset, s_min, s_max, g_Xguard, g_Xfailed, X0_margin, options)

% options
solver_options = spot_sdp_default_options();
solver_options.do_facial_reduction = false;
solver_options.verbose = 0;

bilinear_solve_options.max_iters = 75;
bilinear_solve_options.rank_tol = 1e-5;

use_stored_initial_guess = true;
verify_manual_barrier_function = isfield(options, 'B_manual');
barrier_grow_iters = 5;

% degrees
u_degree = 1;
Nu_degree = 2;
B_degree = 2;
L0_degree = 2;
LF_degree = 2;
NBdot_degree = 2;
s_degree = 2;
Ns_degree = 2;
Ls_guard_degree = 2;
LB_prev_B_degree = 2;
LB_prev_Xfailed_degree = 2;
NB_prev_reset_degree = 2;

% spotless setup
solver = getSolver();
prog = spotsosprog;

if ~isempty(reset)
  % % discrete input
  % ns = length(s_min);
  % [prog, s] = prog.newFreePoly(monomials(x, 0 : s_degree), ns);
  % [prog, Ns_min] = prog.newFreePoly(monomials(x, 0 : Ns_degree));
  % [prog, Ls_guard_min] = prog.newSOSPoly(monomials(x, 0 : Ls_guard_degree));
  % bilinear_sos_constraints{end + 1} = s - s_min + Ns_min * B + Ls_guard_min * g_Xguard(x); % s >= s_min on B(x) = 0, g_Xguard(x) <= 0
  % [prog, Ns_max] = prog.newFreePoly(monomials(x, 0 : Ns_degree));
  % [prog, Ls_guard_max] = prog.newSOSPoly(monomials(x, 0 : Ls_guard_degree));
  % bilinear_sos_constraints{end + 1} = s_max - s + Ns_max * B + Ls_guard_max * g_Xguard(x); % s <= s_max on B(x) = 0, g_Xguard(x) <= 0
  
  % % B_prev(xPrime) <= 0 wherever xPrime = reset(x, s), B(x) <= 0, g_Xguard <= 0
  % [prog, xPrime] = prog.newIndeterminate('y', nstates);
  % [prog, LB_prev_B] = prog.newSOSPoly(monomials(x, 0 : LB_prev_B_degree));
  % [prog, LB_prev_Xfailed] = prog.newSOSPoly(monomials(x, 0 : LB_prev_Xfailed_degree));
  % [prog, NB_prev_reset] = prog.newFreePoly(monomials(x, 0 : NB_prev_reset_degree));
  % bilinear_sos_constraints{end + 1} = -B_prev(xPrime) + NB_prev_reset * reset(x, s) + LB_prev_B * B + LB_prev_Xfailed * g_Xfailed(x);
end

% Indeterminates
[prog, x] = prog.newIndeterminate('x', nstates);

% Inputs
ninputs = length(u_min);
[prog, u] = prog.newFreePoly(monomials(x, 0 : u_degree), ninputs);

% Barrier function
if verify_manual_barrier_function
  B = options.B_manual(x);
  [prog, B_full] = prog.newFreePoly(monomials(x, 0 : B_degree));
  prog = prog.withPolyEqs(B - B_full);
else
  [prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));
end

% Unsafe set constraint
[prog, LF] = prog.newSOSPoly(monomials(x, 0 : LF_degree));
prog = prog.withSOS(B - LF * g_Xfailed(x)); % B >= 0 on g_Xu

% Input limits
bilinear_sos_constraints = cell(0);
[prog, Nu_min] = prog.newFreePoly(monomials(x, 0 : Nu_degree));
bilinear_sos_constraints{end + 1} = u - u_min + Nu_min * B; % u >= u_min on B(x) = 0
[prog, Nu_max] = prog.newFreePoly(monomials(x, 0 : Nu_degree));
bilinear_sos_constraints{end + 1} = u_max - u + Nu_max * B; % u <= u_max on B(x) = 0

% Barrier function derivative constraint
Bdot = diff(B, x) * f(x, u);
[prog, NBdot] = prog.newFreePoly(monomials(x, 0 : NBdot_degree));
epsilon = 0; % TODO
bilinear_sos_constraints{end + 1} = -Bdot + NBdot * B - epsilon; % Bdot <= -epsilon on B(x) = 0

% Initial condition constraint
g_X0 = B_prev(x);
[prog, L0] = prog.newSOSPoly(monomials(x, 0 : L0_degree));
prog = prog.withSOS(-B - L0 * g_X0 - X0_margin); % B <= -X0_margin on g_X0

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
  Bdot_sol = sol.eval(Bdot);
  u_sol = sol.eval(u);
  options.plotfun(B_sol, Bdot_sol, x, u_sol, f);

  % TODO: if success...
  B_fun = makeFunction(B_sol, x);
  u_fun = makeFunction(u_sol, x);
  
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
  prog_base = prog;
  for i = 1 : barrier_grow_iters
    [sol, success, sol_w] = solveBilinear(prog, bilinear_sos_constraints, x, solver, solver_options, bilinear_solve_options);
    disp(['status: ' char(sol.status)]);
    disp(['success: ' num2str(success)])
    fprintf('\n\n');
    
    B_sol = sol.eval(B);
    Bdot_sol = sol.eval(Bdot);
    u_sol = sol.eval(u);
    
    if success
      sol_w_best = sol_w;
      options.plotfun(B_sol, Bdot_sol, x, u_sol, f);
      
      B_fun = makeFunction(B_sol, x);
      u_fun = makeFunction(u_sol, x);
      
      prog = prog_base;
      [prog, L0] = prog.newSOSPoly(monomials(x, 0 : L0_degree));
      prog = prog.withSOS(-B - L0 * -B_sol);
    end

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

function fun = makeFunction(poly, var)
fun = @(input) easySubs(poly, var, input);
end

function ret = easySubs(poly, var, input)
if isnumeric(input)
  ret = full(msubs(poly, var, input));
else
  ret = subs(poly, var, input);
end
end
