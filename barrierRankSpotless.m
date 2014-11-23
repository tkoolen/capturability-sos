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

% Description of initial condition set
g_X0 = 0.5^2 - (x1 - 1.5)^2 - x2^2; % Should be in terms of x1 and x2

% Description of unsafe set
g_Xu = 0.4^2 - (x1+1)^2 - (x2+1)^2; % Should be in terms of x1 and x2

% Multipliers
dL1=2;
[prog, L1] = prog.newSOSPoly(monomials(x, 0 : dL1));

dL2=2;
[prog, L2] = prog.newSOSPoly(monomials(x, 0 : dL2));

dL=2;
[prog, L, cL] = prog.newFreePoly(monomials(x, 0 : dL));

% Barrier function
[prog, B, cB] = prog.newFreePoly(monomials(x, 0 : 4));

% Compute time derivative of B
Bdot = diff(B, x) * f;

constr = -Bdot + L*B; % Bdot < 0 when B = 0
[cconstr, mconstr] = mss_coefficients(constr, x);


w = [cB;cL];
[prog, W] = prog.newSym(length(w));
ww = w*w';

dcconstrdcB = diff(cconstr, cB);
mixed_deriv = sparse(double(diff(dcconstrdcB(:), cL)));
[i, cL_ind] = find(mixed_deriv ~= 0);
coeffs = mixed_deriv(sub2ind(size(mixed_deriv), i, cL_ind));
[cconstr_ind, cB_ind] = ind2sub(size(dcconstrdcB), i);
% at this point, we should have cconstr(cconstr_ind) - cL(cL_ind) .* cB(cB_ind) = 0
cB_w_ind = cB_ind;
cL_w_ind = cL_ind + length(cB);
% at this point, we should have cconstr(cconstr_ind) - ww(sub2ind(size(ww), cB_w_ind, cL_w_ind)) = 0
cconstr_W = cconstr;
cconstr_W(cconstr_ind) = coeffs .* W(sub2ind(size(ww), cB_w_ind, cL_w_ind));
% at this point we should have subs(cconstr_W, mss_s2v(W), mss_s2v(ww)) - cconstr = 0
constr = mconstr * cconstr_W;

M = [W w;w' 1];

% bilinear constraint:
prog = prog.withSOS(constr);
prog = prog.withPSD(M);

% initial condition set constraint
prog = prog.withSOS(-B - L1 * g_X0 - 1e-2);

% unsafe set constraint
prog = prog.withSOS(B - L2 * g_Xu);

sol_w = zeros(size(w));

options = spot_sdp_default_options();
% options.verbose = 1;
for k = 1 : 20 %this loop improves on the standard trace heuristic for rank.
  sol = prog.minimize(trace(W)-2*sol_w'*w,solver,options);
  sol_w = double(sol.eval(w));
%   eig(double(W))%just so we keep track of it.

end
sol_W = double(sol.eval(W));

disp(['solution rank: ' num2str(rank(sol_W, 1e-7))]);

eig_sol_W = eig(sol_W); %if rank 1, we are good, W must be w*w'.
residual = sol_w*sol_w'-sol_W;


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
