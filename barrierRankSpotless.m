function barrierRankSpotless()
solver = getSolver();

checkDependency('spotless');
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

function [prog, sol] = solveBilinear(prog, constr, x, c1, c2, solver, options, max_iters, rank_tol)
% solves SOS program sos subject to additional bilinear SOS constraint
% constr. The constraint constr is bilinear in coefficients c1 and c2, and
% has free variables x.
% From Ibaraki, Tomizuka - Rank Minimization Approach for Solving BMI Problems with Random Search

[cconstr, mconstr] = mss_coefficients(constr, x);
w = [c1; c2];
[prog, W] = prog.newSym(length(w));
ww = w*w';

dcconstrdx1 = diff(cconstr, c1);
mixed_deriv = sparse(double(diff(dcconstrdx1(:), c2)));
[i, c1_ind] = find(mixed_deriv ~= 0);
coeffs = mixed_deriv(sub2ind(size(mixed_deriv), i, c1_ind));
[cconstr_ind, c2_ind] = ind2sub(size(dcconstrdx1), i);

% at this point, we should have cconstr(cconstr_ind) - c2(c2_ind) .* c1(c1_ind) = 0
c1_w_ind = c2_ind;
c2_w_ind = c1_ind + length(c1);
% at this point, we should have cconstr(cconstr_ind) - ww(sub2ind(size(ww), c1_w_ind, c2_w_ind)) = 0
cconstr_W = cconstr;
cconstr_W(cconstr_ind) = coeffs .* W(sub2ind(size(ww), c1_w_ind, c2_w_ind));
% at this point we should have subs(cconstr_W, mss_s2v(W), mss_s2v(ww)) - cconstr = 0
constr = mconstr * cconstr_W;

M = [W w;w' 1];

% bilinear constraint:
prog = prog.withSOS(constr);
prog = prog.withPSD(M);

sol_w = zeros(size(w));

for k = 1 : max_iters %this loop improves on the standard trace heuristic for rank.
  disp(['iteration ' num2str(k)]);
  sol = prog.minimize(trace(W)-2*sol_w'*w,solver,options);
  sol_w = double(sol.eval(w));
  sol_W = double(sol.eval(W));
  W_rank = rank(sol_W, rank_tol);
  if W_rank == 1
    %if rank 1, we are good, W must be w*w'.
    break;
  end
end

disp(['solution rank: ' num2str(W_rank)]);
residual = sol_w * sol_w' - sol_W;
disp(['residual norm: ' num2str(norm(residual))]);
end
