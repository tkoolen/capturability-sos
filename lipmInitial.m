function lipmInitial()
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
[prog, x] = prog.newIndeterminate('x', nstates);
r = x(1);
rd = x(2);

% Dynamics
f = [rd; r];

% Barrier function
dB = 4;
% [prog, B, cB] = prog.newFreePoly(monomials(x, 0 : dB));
% mB = monomialSubs('y', [r + rd; r - rd], 1 : dB);
mB = monomials(x, 1:dB);
[prog, B] = prog.newFreePoly(mB);

% Initial condition constraint
x_star = [0; 0];
A = double(subs(diff(f, x), x, x_star));
[V, D] = eig(A);
unstable_eigenvectors = V(:, diag(D) >= 0);
g_X0 = - (unstable_eigenvectors' * x)' * (unstable_eigenvectors' * x);
% g_X0 = -x' * x;
% g_X0 = - (x - [1; 1])' * (x - [1; 1]);
% g_X0 = -(r + rd)' * (r + rd);

dL1 = 4;
[prog, L1] = prog.newSOSPoly(monomials(x, 0 : dL1));
prog = prog.withSOS(-B - L1 * g_X0);

% Unsafe set constraint
dL2 = 4;
[prog, L2] = prog.newSOSPoly(monomials(x, 0 : dL2));
% rf_dist = 10000;
% g_Xu = r' * r - rf_dist^2;
g_Xu = (r + rd)' * (r + rd) - 0.1;
prog = prog.withSOS(B - L2 * g_Xu - 1);

% Barrier function derivative constraint
Bdot = diff(B, x) * f;
dL = 4;
[prog, L] = prog.newFreePoly(monomials(x, 0 : dL));
constr = -Bdot + L*B; % Bdot < 0 when B = 0

solver_options = spot_sdp_default_options();
solver_options.do_facial_reduction = true;
solver_options.verbose = 0;
options.max_iters = 20;
options.rank_tol = 1e-7;
[prog, sol] = solveBilinear(prog, constr, x, solver, solver_options, options);

B_sol = clean(sol.eval(B), 1e-5);
Bdot_sol = clean(sol.eval(Bdot), 1e-5);

visualize(B_sol, Bdot_sol, g_X0, g_Xu, x, f);

end

function visualize(B_sol, Bdot_sol, g_X0, g_Xu, x, f)
[X,Y] = meshgrid(linspace(-3,3,100),linspace(-3,3,100));
% initial condition set
gs_X0 = msubs(g_X0,x,[X(:),Y(:)]');
% unsafe set
gs_Xu = msubs(g_Xu,x,[X(:),Y(:)]');
gs_B = msubs(B_sol,x,[X(:),Y(:)]');
gs_Bdot = msubs(Bdot_sol,x,[X(:),Y(:)]');

figure();
mesh(X, Y, reshape(double(gs_B),100,100)); xlabel('x1'); ylabel('x2');
title('barrier function');

figure()
mesh(X, Y, reshape(double(gs_Bdot),100,100)); xlabel('x1'); ylabel('x2');
title('barrier function derivative');

figure()
hold on
contour(X,Y,reshape(double(gs_X0),100,100),[0 0],'LineWidth',3); % initial condition set
contour(X,Y,reshape(double(gs_Xu),100,100),[0 0],'r','LineWidth',3) % unsafe set
contour(X,Y,reshape(double(gs_B),100,100),[0 0],'b','LineWidth',3) % boundary function zero level set

% (scaled) vector field
[X,Y] = meshgrid(linspace(-3,3,50),linspace(-3,3,50));
x1dots = reshape(double(msubs(f(1),x,[X(:),Y(:)]')),50,50);
x2dots = reshape(double(msubs(f(2),x,[X(:),Y(:)]')),50,50);
x1dots = 0.1*x1dots./(sqrt(x1dots.^2 + x2dots.^2));
x2dots = 0.1*x2dots./(sqrt(x1dots.^2 + x2dots.^2));
quiver(X,Y,x1dots,x2dots,'AutoScale','off','Color','k');

% Title
title('Barrier functions')
xlabel('x_1');
ylabel('x_2');

legend('X_0','X_u','0 level-set of B(x)');

end
