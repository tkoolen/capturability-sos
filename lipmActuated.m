function lipmActuated()
addpath(fullfile('util'));
cleaner = onCleanup(@() rmpath(fullfile('util')));
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

solver = getSolver();
% solver = @spot_sedumi;
prog = spotsosprog;

% Variables/Indeterminates
nstates = 2;
ninputs = 1;
[prog, x] = prog.newIndeterminate('x', nstates);
[prog, u] = prog.newIndeterminate('u', ninputs);
r = x(1);
rd = x(2);

% input limits
u_max = 0.1;
u_vertices = [-u_max u_max];

% dynamics
f = @(x, u) lipmDynamics(x, u);

% Barrier function
dB = 2;
% [prog, B, cB] = prog.newFreePoly(monomials(x, 0 : dB));
% mB = monomialSubs('y', [r + rd; r - rd], 1 : dB);
mB = monomials(x, 0:dB);
[prog, B] = prog.newFreePoly(mB);

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

dL1 = 3;
[prog, L1] = prog.newSOSPoly(monomials(x, 0 : dL1));
prog = prog.withSOS(-B - L1 * g_X0 - 1e-4); % B <= -1 on g_X0

% Unsafe set constraint
dL2 = 3;
[prog, L2] = prog.newSOSPoly(monomials(x, 0 : dL2));
% rf_dist = 10000;
% g_Xu = r' * r - rf_dist^2;
g_Xu = (r + rd)' * (r + rd) - 10;
prog = prog.withSOS(B - L2 * g_Xu);% - 1); % B >= 1 on g_Xu

% Barrier function derivative constraint
n_u_vertices = size(u_vertices, 2);
Bdot = cell(n_u_vertices, 1);
dB = diff(B, x);
for i = 1 : n_u_vertices
  ui = u_vertices(:, i);
  Bdot{i} = dB * f(x, ui);
end

bilinear_sos_constraints = cell(size(u_vertices, 2), 1);
L_degree = 3;
for j = 1 : size(u_vertices, 2)
  region = msspoly(0);
  for i = 1 : size(u_vertices, 2)
    if i ~= j
      [prog, Lij] = prog.newSOSPoly(monomials(x, 0 : L_degree));
      region = region + Lij * (Bdot{j} - Bdot{i});
    end
  end
  
  [prog, N] = prog.newFreePoly(monomials(x, 0 : L_degree));
  region = region + N * B;
  bilinear_sos_constraints{j} = -Bdot{j} + region;
end

solver_options = spot_sdp_default_options();
solver_options.do_facial_reduction = true;
solver_options.verbose = 0;
options.max_iters = 15;
options.rank_tol = 1e-3;
[prog, sol] = solveBilinear(prog, bilinear_sos_constraints, x, solver, solver_options, options);

B_sol = sol.eval(B);
Bdot_sol = cellfun(@(x) sol.eval(x), Bdot, 'UniformOutput', false);

visualize(B_sol, Bdot_sol, u_vertices, g_X0, g_Xu, x, f);

end


