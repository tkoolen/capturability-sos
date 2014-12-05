function B = capturabilityBarrierRecursive(B_prev, f, nstates, u_min, u_max, reset, s_min, s_max, g_Xguard, g_Xfailed, options)

% parameters
u_degree = 3;
s_degree = 2;
B_degree = 2;
Lfailed_degree = 2;
Nu_degree = 2;
Ns_degree = 2;
Ls_guard_degree = 2;
LB_prev_B_degree = 2;
LB_prev_Xfailed_degree = 2;
NBdot_degree = 2;
NB_prev_reset_degree = 2;

bilinear_sos_constraints = cell(0);
prog = spotsosprog;

% states
[prog, x] = prog.newIndeterminate('x', nstates);

% barrier function
[prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));

% continuous input
nu = length(u_min);
[prog, u] = prog.newFreePoly(monomials(x, 0 : u_degree), nu);
[prog, Nu_min] = prog.newFreePoly(monomials(x, 0 : Nu_degree));
bilinear_sos_constraints{end + 1} = u - u_min + Nu_min * B; % u >= u_min on B(x) = 0
[prog, Nu_max] = prog.newFreePoly(monomials(x, 0 : Nu_degree));
bilinear_sos_constraints{end + 1} = u_max - u + Nu_max * B; % u <= u_max on B(x) = 0

% discrete input
ns = length(s_min);
[prog, s] = prog.newFreePoly(monomials(x, 0 : s_degree), ns);
[prog, Ns_min] = prog.newFreePoly(monomials(x, 0 : Ns_degree));
[prog, Ls_guard_min] = prog.newSOSPoly(monomials(x, 0 : Ls_guard_degree));
bilinear_sos_constraints{end + 1} = s - s_min + Ns_min * B + Ls_guard_min * g_Xguard(x); % s >= s_min on B(x) = 0, g_Xguard(x) <= 0
[prog, Ns_max] = prog.newFreePoly(monomials(x, 0 : Ns_degree));
[prog, Ls_guard_max] = prog.newSOSPoly(monomials(x, 0 : Ls_guard_degree));
bilinear_sos_constraints{end + 1} = s_max - s + Ns_max * B + Ls_guard_max * g_Xguard(x); % s <= s_max on B(x) = 0, g_Xguard(x) <= 0

% Unsafe set constraint
[prog, Lfailed] = prog.newSOSPoly(monomials(x, 0 : Lfailed_degree));
prog = prog.withSOS(B - Lfailed * g_Xfailed(x)); % B >= 0 on Xfailed

% B_prev(xPrime) <= 0 wherever xPrime = reset(x, s), B(x) <= 0, g_Xguard <= 0
[prog, xPrime] = prog.newIndeterminate('y', nstates);
[prog, LB_prev_B] = prog.newSOSPoly(monomials(x, 0 : LB_prev_B_degree));
[prog, LB_prev_Xfailed] = prog.newSOSPoly(monomials(x, 0 : LB_prev_Xfailed_degree));
[prog, NB_prev_reset] = prog.newFreePoly(monomials(x, 0 : NB_prev_reset_degree));
bilinear_sos_constraints{end + 1} = -B_prev(xPrime) + NB_prev_reset * reset(x, s) + LB_prev_B * B + LB_prev_Xfailed * g_Xfailed(x);

% Barrier function derivative constraint
Bdot = diff(B, x) * f(x, u); % TODO
[prog, NBdot] = prog.newFreePoly(monomials(x, 0 : NBdot_degree));
epsilon = 1e-10;
bilinear_sos_constraints{end + 1} = -Bdot + NBdot * B - epsilon; % Bdot <= -epsilon on B(x) = 0

[sol, success, sol_w] = solveBilinear(prog, bilinear_sos_constraints, [x; xPrime], options.solver, options.solver_options, options.bilinear_solve_options);
disp(['status: ' char(sol.status)]);
disp(['success: ' num2str(success)])
fprintf('\n\n');

end