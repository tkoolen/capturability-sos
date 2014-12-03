function [prog, sol] = solveBilinear(prog, constraints, x, solver, solver_options, options)
% solves SOS program prog subject to additional bilinear SOS constraint
% constr. The constraint constr is bilinear in coefficients c1 and c2, and
% has free variables x.
% From Ibaraki, Tomizuka - Rank Minimization Approach for Solving BMI Problems with Random Search

if ~isfield(options, 'max_iters')
  options.max_iters = 20;
end
if ~isfield(options, 'rank_tol')
  options.rank_tol = 1e-5;
end
if ~isfield(options, 'objective_type')
  options.objective_type = 'sum';
end
if ~isfield(options, 'use_small_psd_constraints')
  options.use_small_psd_constraints = true;
end
 
if ~iscell(constraints)
  constraints = {constraints};
end

n_constraints = length(constraints);
W = cell(n_constraints, 1);
w = cell(n_constraints, 1);
M = cell(n_constraints, 1);
sol_w = cell(n_constraints, 1);

for j = 1 : n_constraints
  constraint = constraints{j};
  
  % decompose constraint
  [vars, degrees, monomials] = decomp(constraint, x);
  coefficients = recomp(vars, degrees, speye(size(degrees, 1)));
  
  total_degrees = sum(degrees, 2);
  nonlinear_monomial_indices = find(total_degrees > 1);
  n_nonlinear_monomial_indices = length(nonlinear_monomial_indices);
  
  if options.use_small_psd_constraints % Create symmetric matrix + PSD constraint for each bilinear term in each constraint
    W{j} = cell(n_nonlinear_monomial_indices, 1);
    w{j} = cell(n_nonlinear_monomial_indices, 1);
    M{j} = cell(n_nonlinear_monomial_indices, 1);
    sol_w{j} = cell(n_nonlinear_monomial_indices, 1);
   
    coefficients_linear = coefficients;
    for i = 1 : n_nonlinear_monomial_indices
      coeff_idx = nonlinear_monomial_indices(i);
      [prog, coefficients_linear(coeff_idx), w{j, i}, W{j, i}, M{j, i}] = replaceBilinearTermsWithNewVariables(prog, coefficients(coeff_idx));
      sol_w{j, i} = zeros(size(w{j, i}));
    end
  else % create one symmetric matrix + PSD constraint for each constraint
    [prog, coefficients_linear, w{j}, W{j}, M{j}] = replaceBilinearTermsWithNewVariables(prog, coefficients);
    sol_w{j} = zeros(size(w{j}));
  end
  
  % recompute constraint in terms of bilinear variables
  constr_linear = monomials * coefficients_linear;
  
  prog = prog.withSOS(constr_linear);
end

% iteratively solve
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
  sol_w = cellfun(@(x) double(sol.eval(x)), w, 'UniformOutput', false);
  M_ranks = cellfun(@(x) rank(double(sol.eval(x)), options.rank_tol), M);
  
  disp(['objective value: ' num2str(double(sol.eval(f)))]);
  disp(['max rank: ' num2str(max(max(M_ranks)))]);
  disp(['number of bilinear variable matrices with rank > 1: ' num2str(sum(M_ranks(:) > 1))]);
  fprintf('\n');
  
  if max(max(M_ranks)) == 1
    break;
  end
end


% residual = sol_w * sol_w' - sol_W;
% disp(['residual norm: ' num2str(norm(residual))]);

end

function f = sumObjective(W, w, sol_w)
f = zeros(1, 1, 'msspoly');
for i = 1 : numel(W)
  f = f + trace(W{i}) - 2 * sol_w{i}' * w{i};
end
end
