function lipm()
solver = getSolver();

checkDependency('spotless');
prog = spotsosprog;

nstates = 2;
ninputs = 1;

[prog, x] = prog.newIndeterminate('x', nstates);
r = x(1);
rd = x(2);

% dynamics
f = @(x, u) lipmDynamics(x, u);

% input limits
u_max = 0.1;
u_vertices = [-u_max u_max];

% barrier function
B_degree = 4;
[prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));

% failed states
rf_dist = 10000;
g_Xf = r' * r - rf_dist^2;

% fixed point
x_fixed = [0; 0];
u_fixed = 0;

% initial states
g_X0 = -(x - x_fixed)' * (x - x_fixed);
% g_X0 = -(x(1) + x(2))^2;

L_degree = 2;

% [prog, L0] = prog.newSOSPoly(monomials(x, 0 : L_degree));
% prog = prog.withSOS(-B - L0 * g_X0 - 1e-4);
prog = prog.withSOS(subs(-B - 1e-4, x, x_fixed));

[prog, Lf] = prog.newSOSPoly(monomials(x, 0 : L_degree));
prog = prog.withSOS(B - Lf * g_Xf);

Bdot_fixed = diff(B, x) * f(x, u_fixed);
prog = prog.withSOS(-Bdot_fixed);

% n_u_vertices = size(u_vertices, 2);
% Bdot = cell(n_u_vertices, 1);
% dB = diff(B, x);
% for i = 1 : n_u_vertices
%   ui = u_vertices(:, i);
%   Bdot{i} = dB * ui;
% end
% 
% for j = 1 : size(u_vertices, 2)
%   region = msspoly(0);
%   for i = 1 : size(u_vertices, 2)
%     if i ~= j
%       [prog, Lij] = prog.newSOSPoly(monomials(x, 0 : L_degree));
%       
%       % can't do this because both Lij and Bdot{j} - Bdot{i} are free
%       % polynomials, so product is not linear in coefficients:
% %       region = region + Lij * (Bdot{j} - Bdot{i});
%     end
%   end
%   
%   prog = prog.withSOS(-Bdot{j} + region);
% end

options = spot_sdp_default_options();
options.verbose = 1;
sol = prog.minimize(0,solver,options);

% Check is SOS program ran correctly
if ~sol.isPrimalFeasible
  error('The SOS problem is not feasible');
end

B_sol = clean(sol.eval(B),1e-5);
% disp('Barrier function:');
% B_sol

% Plot stuff
figure
[X,Xd] = meshgrid(linspace(-3,3,100),linspace(-3,3,100));
% initial condition set
gs_X0 = msubs(g_X0,x,[X(:),Xd(:)]');
contour(X,Xd,reshape(double(gs_X0),100,100),[0 0],'g','LineWidth',3); % initial condition set
hold on
% unsafe set
gs_Xu = msubs(g_Xf,x,[X(:),Xd(:)]');
contour(X,Xd,reshape(double(gs_Xu),100,100),[0 0],'r','LineWidth',3) % unsafe set
% 0-level set of B
gs_B = msubs(sol.eval(B),x,[X(:),Xd(:)]');
contour(X,Xd,reshape(double(gs_B),100,100),[0 0],'b','LineWidth',3) % barrier zero level set

% (scaled) vector field
% [X,Xd] = meshgrid(linspace(-3,3,50),linspace(-3,3,50));
% x1dots = reshape(double(msubs(f(1),x,[X(:),Xd(:)]')),50,50);
% x2dots = reshape(double(msubs(f(2),x,[X(:),Xd(:)]')),50,50);
% x1dots = 0.1*x1dots./(sqrt(x1dots.^2 + x2dots.^2));
% x2dots = 0.1*x2dots./(sqrt(x1dots.^2 + x2dots.^2));
% quiver(X,Xd,x1dots,x2dots,'AutoScale','off','Color','k');

% Title
title('Barrier functions')
xlabel('x_1');
ylabel('x_2');
legend('X_0','X_u','0 level-set of B(x)');

end
