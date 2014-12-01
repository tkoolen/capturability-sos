function [prog, sol] = solveBilinear(prog, constr, x, solver, options, max_iters, rank_tol)
% solves SOS program prog subject to additional bilinear SOS constraint
% constr. The constraint constr is bilinear in coefficients c1 and c2, and
% has free variables x.
% From Ibaraki, Tomizuka - Rank Minimization Approach for Solving BMI Problems with Random Search

[vars, degrees, monomials] = decomp(constr, x);
n_vars = length(vars);
total_degrees = sum(degrees, 2);
if any(total_degrees > 2)
  error('constraint is not bilinear')
end

bilinear_indices = find(total_degrees == 2);
degrees_bilinear = degrees(bilinear_indices, :);
n_bilinear_vars = size(degrees_bilinear, 1);
bilinear_vars = zeros(n_bilinear_vars, 1, 'msspoly');

% create symmetric matrix + PSD constraint for each monomial that is bilinear
W = cell(n_bilinear_vars, 1);
w = cell(n_bilinear_vars, 1);
sol_w = cell(n_bilinear_vars, 1);

[prog, Wvec] = prog.newFree(spotprog.psdDimToNo(2), n_bilinear_vars);

for i = 1 : n_bilinear_vars
  w{i} = vars(logical(degrees_bilinear(i, :)));
  W{i} = mss_v2s(Wvec(:, i));
  M = [W{i} w{i}; w{i}' 1];
  prog = prog.withPSD(M);
  bilinear_vars(i) = W{i}(1, 2);
  sol_w{i} = zeros(size(w{i}));
end

% change degrees so that we use the bilinear variables instead of products
% of the original variables
degrees(bilinear_indices, :) = 0;
[i,j,v] = find(degrees);
i = [i; bilinear_indices];
j = [j; n_vars + (1 : n_bilinear_vars)'];
v = [v; ones(size(bilinear_indices))];
degrees = sparse(i, j, v, size(degrees, 1), n_vars + n_bilinear_vars);

% recompose coefficients after getting rid of bilinear monomials and adding
% the bilinear variables instead
coefficients_linear = recomp([vars; bilinear_vars], degrees, speye(size(degrees, 1)));

% recompute constraint in terms of bilinear variables
constr_linear = monomials * coefficients_linear;

constr_back = constr_linear;
for i = 1 : n_bilinear_vars
  constr_back = subs(constr_back, bilinear_vars(i), prod(w{i}));
end
valuecheck(0, double(constr_back - constr));

prog = prog.withSOS(constr_linear);

for k = 1 : max_iters %this loop improves on the standard trace heuristic for rank.
  disp(['iteration ' num2str(k)]);
  f = objective(W, w, sol_w);
  sol = prog.minimize(f, solver,options);
  sol_w = getSolW(sol, w);
  if double(sol.eval(f)) < rank_tol
    break;
  end
end

disp(['objective value: ' num2str(double(sol.eval(f)))]);
% residual = sol_w * sol_w' - sol_W;
% disp(['residual norm: ' num2str(norm(residual))]);

end

function f = objective(W, w, sol_w)
f = zeros(1, 1, 'msspoly');
for i = 1 : length(W)
  f = f + trace(W{i}) - 2 * sol_w{i}' * w{i};
end
end

function sol_w = getSolW(sol, w)
sol_w = cell(size(w));
for i = 1 : length(w)
  sol_w{i} = double(sol.eval(w{i}));
end
end