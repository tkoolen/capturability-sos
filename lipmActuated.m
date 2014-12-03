function lipmActuated()
% path setup
addpath(fullfile('util'));
cleaner = onCleanup(@() rmpath(fullfile('util')));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

solver = getSolver();
solver_options = spot_sdp_default_options();
solver_options.do_facial_reduction = false;
solver_options.verbose = 0;
bilinear_solve_options.max_iters = 50;
bilinear_solve_options.rank_tol = 1e-3;
verify_manual_barrier_function = false;
L0_degree = 3;
Lu_degree = 3;
Lij_degree = 3;
N_degree = 3;

prog = spotsosprog;

% Variables/Indeterminates
nstates = 2;
ninputs = 1;
[prog, x] = prog.newIndeterminate('x', nstates);
[prog, u] = prog.newIndeterminate('u', ninputs);
r = x(1);
rd = x(2);

% input limits
u_max = 1;
u_vertices = [-u_max u_max];

% dynamics
f = @(x, u) lipmDynamics(x, u);

% Barrier function
if verify_manual_barrier_function
  B = (r + rd)^2 / u_max^2 - 1;
else
  dB = 3;
  [prog, B] = prog.newFreePoly(monomials(x, 0 : dB));
end

% Initial condition constraint
x_star = [0; 0];
u_star = 0;
A = double(subs(diff(f(x, u), x), [x; u], [x_star; u_star]));
[V, D] = eig(A);
unstable_eigenvectors = V(:, diag(D) >= 0);
g_X0 = - (unstable_eigenvectors' * x)' * (unstable_eigenvectors' * x);
% g_X0 = - (x - x_star)' * (x - x_star);
% g_X0 = -x' * x;
% g_X0 = - (x - [1; 1])' * (x - [1; 1]);
% g_X0 = -(r + rd)' * (r + rd);

[prog, L0] = prog.newSOSPoly(monomials(x, 0 : L0_degree));
prog = prog.withSOS(-B - L0 * g_X0 - 1); % B <= -1 on g_X0

% Unsafe set constraint
[prog, Lu] = prog.newSOSPoly(monomials(x, 0 : Lu_degree));
% rf_dist = 10000;
% g_Xu = r' * r - rf_dist^2;
x_ic_dist = 0.5;
g_Xu = (r + rd)' * (r + rd) - x_ic_dist^2;
prog = prog.withSOS(B - Lu * g_Xu);% - 1); % B >= 1 on g_Xu

% Barrier function derivative constraint
n_u_vertices = size(u_vertices, 2);
Bdot = cell(n_u_vertices, 1);
dB = diff(B, x);
for i = 1 : n_u_vertices
  ui = u_vertices(:, i);
  Bdot{i} = dB * f(x, ui);
end

bilinear_sos_constraints = cell(size(u_vertices, 2), 1);
for j = 1 : size(u_vertices, 2)
  region = msspoly(0);
  for i = 1 : size(u_vertices, 2)
    if i ~= j
      [prog, Lij] = prog.newSOSPoly(monomials(x, 0 : Lij_degree));
      region = region + Lij * (Bdot{j} - Bdot{i});
    end
  end
  
  [prog, N] = prog.newFreePoly(monomials(x, 0 : N_degree));
  region = region + N * B;
  
  if verify_manual_barrier_function
    prog = prog.withSOS(-Bdot{j} + region);
  else
    bilinear_sos_constraints{j} = -Bdot{j} + region;
  end
end

if verify_manual_barrier_function
  sol = prog.minimize(0, solver, solver_options);
  disp(['status: ' char(sol.status)]);
else
  [~, sol] = solveBilinear(prog, bilinear_sos_constraints, x, solver, solver_options, bilinear_solve_options);
end

B_sol = sol.eval(B);
Bdot_sol = cellfun(@(x) sol.eval(x), Bdot, 'UniformOutput', false);

visualize(B_sol, Bdot_sol, u_vertices, g_X0, g_Xu, x, f);

end


