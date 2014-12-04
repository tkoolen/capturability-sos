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
B_degree = 2;
L0_degree = 2;
Lu_degree = 2;
Lij_degree = 2;
N_degree = 2;

verify_manual_barrier_function = false;
use_fixed_control_law = false;

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
  B = (r + rd)^2 / (u_max)^2 - 1;
  [prog, B_full] = prog.newFreePoly(monomials(x, 0 : B_degree));
  prog = prog.withPolyEqs(B - B_full);
else
  [prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));
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
x_ic_dist = 3;
g_Xu = (r + rd)' * (r + rd) - x_ic_dist^2;
prog = prog.withSOS(B - Lu * g_Xu);% - 1); % B >= 1 on g_Xu

% Barrier function derivative constraint
n_u_vertices = size(u_vertices, 2);
Bdot = cell(n_u_vertices, 1);
dB = diff(B, x);
for j = 1 : n_u_vertices
  ui = u_vertices(:, j);
  Bdot{j} = dB * f(x, ui);
end

Ni = cell(n_u_vertices, 1);
Lij = cell(n_u_vertices, n_u_vertices);
if use_fixed_control_law
  A = [1 u_max u_max^2 u_max^3;
    1 -u_max u_max^2 -u_max^3;
    0 1 2 * u_max 3 * u_max^2;
    0 1 -2 * u_max 3 * u_max^2];
  b = [u_max; -u_max; 0; 0];
  h_coeffs = A \ b;
  h_monomials = monomialSubs('z', r + rd, 0:3);
  h = h_coeffs' * h_monomials;
%   prog = prog.withSOS(-dB * f(x, h));
  
  [prog, N] = prog.newFreePoly(monomials(x, 0 : N_degree));
  region = N * B;
  bilinear_sos_constraints = -dB * f(x, h) + region;
else
  bilinear_sos_constraints = cell(size(u_vertices, 2), 1);
  for i = 1 : size(u_vertices, 2)
    region = msspoly(0);
    for j = 1 : size(u_vertices, 2)
      if j ~= i
        [prog, Lij{i, j}] = prog.newSOSPoly(monomials(x, 0 : Lij_degree));
        region = region - Lij{i, j} * (Bdot{j} - Bdot{i});
      end
    end
  
    [prog, Ni{i}] = prog.newFreePoly(monomials(x, 0 : N_degree));
    region = region + Ni{i} * B;
  
    if verify_manual_barrier_function
      prog = prog.withSOS(-Bdot{i} + region);
    else
      bilinear_sos_constraints{i} = -Bdot{i} + region;
    end
  end
end
  
if verify_manual_barrier_function
  sol = prog.minimize(0, solver, solver_options);
else
  [~, sol] = solveBilinear(prog, bilinear_sos_constraints, x, solver, solver_options, bilinear_solve_options);
end

disp(['status: ' char(sol.status)]);

B_sol = sol.eval(B);
Bdot_sol = cellfun(@(x) sol.eval(x), Bdot, 'UniformOutput', false);

visualize(B_sol, Bdot_sol, u_vertices, g_X0, g_Xu, x, f);

end


