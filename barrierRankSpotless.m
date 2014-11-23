function barrierRankSpotless()
addpath(fullfile('util'));
cleaner = onCleanup(@() rmpath(fullfile('util')));
checkDependency('spotless');
solver = getSolver();

prog = spotsosprog;

% Variables/Indeterminates
nstates = 2;
[prog, x] = prog.newIndeterminate('x', nstates);
x1 = x(1);
x2 = x(2);

% Dynamics
f = [x2; -x1 + x1^3 - x2];

% Barrier function
dB = 4;
[prog, B, cB] = prog.newFreePoly(monomials(x, 0 : dB));

% Initial condition constraint
g_X0 = 0.5^2 - (x1 - 1.5)^2 - x2^2; % Should be in terms of x1 and x2
dL1=2;
[prog, L1] = prog.newSOSPoly(monomials(x, 0 : dL1));
prog = prog.withSOS(-B - L1 * g_X0 - 1e-2);

% Unsafe set constraint
dL2=2;
[prog, L2] = prog.newSOSPoly(monomials(x, 0 : dL2));
g_Xu = 0.4^2 - (x1+1)^2 - (x2+1)^2; % Should be in terms of x1 and x2
prog = prog.withSOS(B - L2 * g_Xu);

% Barrier function derivative constraint
Bdot = diff(B, x) * f;
dL=2;
[prog, L, cL] = prog.newFreePoly(monomials(x, 0 : dL));
constr = -Bdot + L*B; % Bdot < 0 when B = 0

options = spot_sdp_default_options();
% options.verbose = 1;
max_iters = 20;
rank_tol = 1e-7;
[prog, sol] = solveBilinear(prog, constr, x, cB, cL, solver, options, max_iters, rank_tol);

% Print out B after zeroing out very small coefficient
xprint = msspoly('x',2); % print in terms of x1 and x2
disp(' ');
disp('Barrier function:');
B_sol = clean(subs(sol.eval(B),x,xprint),1e-5)

% Plot stuff
figure
[X,Y] = meshgrid(linspace(-3,3,100),linspace(-3,3,100));
% initial condition set
gs_X0 = msubs(g_X0,x,[X(:),Y(:)]');
contour(X,Y,reshape(double(gs_X0),100,100),[0 0],'LineWidth',3); % initial condition set
hold on
% unsafe set
gs_Xu = msubs(g_Xu,x,[X(:),Y(:)]');
contour(X,Y,reshape(double(gs_Xu),100,100),[0 0],'r','LineWidth',3) % unsafe set
% 0-level set of B
gs_B = msubs(sol.eval(B),x,[X(:),Y(:)]');
contour(X,Y,reshape(double(gs_B),100,100),[0 0],'b','LineWidth',3) % unsafe set

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
