function [sol, success, w, sol_w] = solveBilinear(prog, constraints, x, solver, solver_options, options)
% solves SOS program prog subject to additional bilinear SOS constraint
% constr. The constraint constr is bilinear in coefficients c1 and c2, and
% has free variables x.
% From Ibaraki, Tomizuka - Rank Minimization Approach for Solving BMI Problems with Random Search

import bilinearSOS.*

if ~isfield(options, 'max_iters')
  options.max_iters = 20;
end
if ~isfield(options, 'rank_tol')
  options.rank_tol = 1e-5;
end
if ~isfield(options, 'tol')
  options.tol = 1e-4;
end
if ~isfield(options, 'objective_type')
  options.objective_type = 'sum';
end
if ~isfield(options, 'psd_constraint_size')
  options.psd_constraint_size = 3;
end

if ~iscell(constraints)
  constraints = {constraints};
end

n_constraints = length(constraints);
coefficients = cell(n_constraints, 1);
monomials = cell(n_constraints, 1);

vars_original = prog.variables;

for j = 1 : n_constraints
  % find coefficients of constraint
  constraint = constraints{j};
  [vars, degrees, monomials{j}] = decomp(constraint, x);
  coefficients{j} = recomp(vars, degrees, speye(size(degrees, 1)));
end

if options.psd_constraint_size == 3 % create one symmetric matrix + PSD constraint for all constraints
  [prog, coefficients_linear_cat, w{1}, W{1}, M{1}] = replaceBilinearTermsWithNewVariables(prog, vertcat(coefficients{:}));
  rows = cellfun(@(x) size(x, 1), coefficients);
  coefficients_linear = mat2cell(coefficients_linear_cat, rows);
else
  W = cell(n_constraints, 1);
  w = cell(n_constraints, 1);
  M = cell(n_constraints, 1);
  coefficients_linear = cell(n_constraints, 1);
  for j = 1 : n_constraints
    if options.psd_constraint_size == 2 % create one symmetric matrix + PSD constraint for each constraint
      [prog, coefficients_linear{j}, w{j}, W{j}, M{j}] = replaceBilinearTermsWithNewVariables(prog, coefficients{j});
    elseif options.psd_constraint_size == 1 % Create symmetric matrix + PSD constraint for each bilinear term in each constraint
      total_degrees = sum(degrees, 2);
      nonlinear_monomial_indices = find(total_degrees > 1);
      n_nonlinear_monomial_indices = length(nonlinear_monomial_indices);
      W{j} = cell(n_nonlinear_monomial_indices, n_nonlinear_monomial_indices);
      w{j} = cell(n_nonlinear_monomial_indices, n_nonlinear_monomial_indices);
      M{j} = cell(n_nonlinear_monomial_indices, n_nonlinear_monomial_indices);
      coefficients_linear{j} = coefficients{j};
      for i = 1 : n_nonlinear_monomial_indices
        coeff_idx = nonlinear_monomial_indices(i);
        [prog, coefficients_linear{j}(coeff_idx), w{j, i}, W{j, i}, M{j, i}] = replaceBilinearTermsWithNewVariables(prog, coefficients{j}(coeff_idx));
      end
    else
      error('psd_constraint_size not recognized');
    end
  end
end

% initial guess
sol_w = cell(size(w));
if isfield(options, 'initial_guess') % initial_guess has values for all vars in prog before adding bilinear vars, so:
  for i = 1 : size(w, 1)
    for j = 1 : size(w, 2)
      bilinear_var_indices = match(vars_original, w{i, j});
      sol_w{i, j} = options.initial_guess(bilinear_var_indices);
    end
  end
else
  for i = 1 : size(w, 1)
    for j = 1 : size(w, 2)
      sol_w{i, j} = 0.1 * randn(size(w{i, j}));
    end
  end
end

% add SOS constraints using linear coefficients
for j = 1 : length(constraints)
  constr_linear = monomials{j} * coefficients_linear{j};
  prog = prog.withSOS(constr_linear);
end

% iteratively solve
t = inf;
success = false;
for k = 1 : options.max_iters
  disp(['iteration ' num2str(k)]);
  if strcmp(options.objective_type, 'sum')
    f = sumObjective(W, w, sol_w);
    prog_k = prog;
  elseif strcmp(options.objective_type, 'max')
    prog_k = prog;
    [prog_k, gamma] = prog_k.newFree(1);
    for i = 1 : numel(W)
      prog_k = prog_k.withPos(gamma - (trace(W{i}) - 2 * sol_w{i}' * w{i}));
    end
    f = gamma;
  else
    error('objective type not recognized');
  end
  
  sol = prog_k.minimize(f, solver, solver_options);
  
  sol_w_new = cellfun(@(x) double(sol.eval(x)), w, 'UniformOutput', false);
  sol_W = cellfun(@(x) double(sol.eval(x)), W, 'UniformOutput', false);
  M_ranks = cellfun(@(x) rank(double(sol.eval(x)), options.rank_tol), M);
  t_new = cellfun(@(W, w, w_new) trace(W) - 2 * w' * w_new + w' * w, sol_W, sol_w, sol_w_new);
  delta_t = t - t_new;
  sol_w = sol_w_new;
  t = t_new;

  disp(['max rank: ' num2str(max(max(M_ranks)))]);
  disp(['number of bilinear variable matrices with rank > 1: ' num2str(sum(M_ranks(:) > 1))]);
  disp(['max(t): ' num2str(max(t))]);
  disp(['max(delta t): ' num2str(max(delta_t))]);
 
  if all(abs(delta_t) < options.tol)
    if t < 10 * options.tol % TODO: make separate parameter
      success = true;
    end
    break;
  end

  fprintf('\n');
end

end

function f = sumObjective(W, w, sol_w)
f = zeros(1, 1, 'msspoly');
for i = 1 : numel(W)
  f = f + trace(W{i}) - 2 * sol_w{i}' * w{i};
end
end

function Sigma = computeCovariance(w, W, alpha)
Sigma = alpha * (W - w * w');
[V, D] = eig(Sigma);
d = diag(D);
d(d <= 0) = 1e-10;
Sigma = (V * diag(d)) / V;
end
