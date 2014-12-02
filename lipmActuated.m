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
dB = 3;
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
prog = prog.withSOS(-B - L1 * g_X0 - 1e-3); % B <= -1 on g_X0

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
options.max_iters = 40;
options.rank_tol = 1e-3;
[prog, sol] = solveBilinear(prog, bilinear_sos_constraints, x, solver, solver_options, options);

B_sol = sol.eval(B);
Bdot_sol = cellfun(@(x) sol.eval(x), Bdot, 'UniformOutput', false);

visualize(B_sol, Bdot_sol, g_X0, g_Xu, x, f);

end

function visualize(B_sol, Bdot_sol, g_X0, g_Xu, x, f)
nmesh = 300;

[X,Y] = meshgrid(linspace(-3,3,nmesh),linspace(-3,3,nmesh));
% initial condition set
gs_X0 = double(msubs(g_X0,x,[X(:),Y(:)]'));
% unsafe set
gs_Xu = double(msubs(g_Xu,x,[X(:),Y(:)]'));
gs_B = double(msubs(B_sol,x,[X(:),Y(:)]'));
gs_Bdot = cellfun(@(z) double(msubs(z,x,[X(:),Y(:)]')), Bdot_sol, 'UniformOutput', false);
gs_Bdotmin = min(vertcat(gs_Bdot{:}), [], 1);

figure();
mesh(X, Y, reshape(double(gs_B),nmesh,nmesh)); xlabel('x1'); ylabel('x2');
title('barrier function');

figure()
hold on;
% for i = 1 : numel(gs_Bdot)
%   mesh(X, Y, reshape(double(gs_Bdot{i}),100,100)); xlabel('x1'); ylabel('x2');
% end
mesh(X, Y, reshape(gs_Bdot{2} - gs_Bdot{1}, nmesh, nmesh)); xlabel('x1'); ylabel('x2');
title('barrier function derivative difference');

figure();
mesh(X, Y, reshape(gs_Bdotmin, nmesh, nmesh)); xlabel('x1'); ylabel('x2');
title('barrier function min derivative');


figure()
hold on
contour(X,Y,reshape(gs_X0,nmesh,nmesh),[0 0],'LineWidth',3); % initial condition set
contour(X,Y,reshape(gs_Xu,nmesh,nmesh),[0 0],'r','LineWidth',3) % unsafe set
contour(X,Y,reshape(gs_B,nmesh,nmesh),[0 0],'b','LineWidth',3) % boundary function zero level set
contour(X,Y,reshape(gs_Bdot{2} - gs_Bdot{1},nmesh,nmesh),[0 0],'m','LineWidth',3) % boundary function derivative difference zero level set


% (scaled) vector field
% [X,Y] = meshgrid(linspace(-3,3,50),linspace(-3,3,50));
% x1dots = reshape(double(msubs(f(1),x,[X(:),Y(:)]')),50,50);
% x2dots = reshape(double(msubs(f(2),x,[X(:),Y(:)]')),50,50);
% x1dots = 0.1*x1dots./(sqrt(x1dots.^2 + x2dots.^2));
% x2dots = 0.1*x2dots./(sqrt(x1dots.^2 + x2dots.^2));
% quiver(X,Y,x1dots,x2dots,'AutoScale','off','Color','k');

% Title
title('Barrier functions')
xlabel('x_1');
ylabel('x_2');

legend('X_0','X_u','0 level-set of B(x)');

end
