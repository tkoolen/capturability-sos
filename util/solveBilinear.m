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
  options.use_small_psd_constraints = false;
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
  n_vars = length(vars);
  total_degrees = sum(degrees, 2);
  if any(total_degrees > 2)
    error('constraint is not bilinear')
  end
  
  bilinear_indices = find(total_degrees == 2);
  degrees_bilinear = degrees(bilinear_indices, :);
  n_bilinear_vars = size(degrees_bilinear, 1);
  bilinear_vars = zeros(n_bilinear_vars, 1, 'msspoly');
  
  if options.use_small_psd_constraints
    % Create symmetric matrix + PSD constraint for each bilinear term.
    W{j} = cell(n_bilinear_vars, 1);
    w{j} = cell(n_bilinear_vars, 1);
    M{j} = cell(n_bilinear_vars, 1);
    sol_w{j} = cell(n_bilinear_vars, 1);
    
    % Create all variables for the symmetric W matrices at once instead of in
    % separate newSymmetric calls inside the loop for efficiency
    [prog, Wvec] = prog.newFree(spotprog.psdDimToNo(2), n_bilinear_vars);
    
    for i = 1 : n_bilinear_vars
      % reshape Wvec variables into a symmetric matrix
      w{j, i} = vars(logical(degrees_bilinear(i, :)));
      W{j, i} = mss_v2s(Wvec(:, i));
      
      % set up the PSD constraint
      M{j, i} = [W{j, i} w{j, i}; w{j, i}' 1];
      prog = prog.withPSD(M{j, i});
      
      % store the 'bilinear variable' that replaces prod(w{i}) in the
      % constraint
      bilinear_vars(i) = W{j, i}(1, 2);
      
      % initial guess for solution
      sol_w{j, i} = zeros(size(w{j, i}));
    end
    
    % change degrees so that we use the bilinear variables instead of products
    % of the original variables
    degrees(bilinear_indices, :) = 0;
    [row, col, value] = find(degrees);
    row = [row; bilinear_indices]; %#ok<AGROW>
    col = [col; n_vars + (1 : n_bilinear_vars)']; %#ok<AGROW>
    value = [value; ones(size(bilinear_indices))]; %#ok<AGROW>
    degrees = sparse(row, col, value, size(degrees, 1), n_vars + n_bilinear_vars);
    
    % recompose coefficients after getting rid of bilinear monomials and adding
    % the bilinear variables instead
    coefficients_linear = recomp([vars; bilinear_vars], degrees, speye(size(degrees, 1)));
  else
    % figure out in which variables the constraint is bilinear
    w{j} = vars(any(degrees_bilinear > 0, 1));
    
    % initial guess for solution
    sol_w{j} = zeros(size(w{j}));
    
    % create symmetric W matrix
    n_Wvec = spotprog.psdDimToNo(length(w{j}));
    [prog, Wvec] = prog.newFree(n_Wvec, 1);
    W{j} = mss_v2s(Wvec);
    
    % create matrix that can be used to look up indices of elements of W in
    % Wvec
    Wvec_index_lookup = mss_v2s(1 : n_Wvec);
    
    % set up the PSD constraint
    M{j} = [W{j} w{j}; w{j}' 1];
    prog = prog.withPSD(M{j});
    
    % find indices into Wvec corresponding to bilinear variables
    Wvec_indices = zeros(size(degrees_bilinear, 1), 1);
    for i = 1 : size(degrees_bilinear, 1)
      var_indices = find(degrees_bilinear(i, :));
      Wvec_indices(i) = Wvec_index_lookup(var_indices(1), var_indices(2));
    end
    
    % change degrees so that we use the bilinear variables instead of products
    % of the original variables
    degrees(bilinear_indices, :) = 0;
    [row, col, value] = find(degrees);
    row = [row; bilinear_indices]; %#ok<AGROW>
    col = [col; n_vars + Wvec_indices]; %#ok<AGROW>
    value = [value; ones(size(Wvec_indices))]; %#ok<AGROW>
    degrees = sparse(row, col, value, size(degrees, 1), n_vars + length(Wvec));
    
    % recompose coefficients after getting rid of bilinear monomials and adding
    % the bilinear variables instead
    coefficients_linear = recomp([vars; Wvec], degrees, speye(size(degrees, 1)));
  end
  
  % recompute constraint in terms of bilinear variables
  constr_linear = monomials * coefficients_linear;
  
  % check
%   constr_back = constr_linear;
%   for i = 1 : n_bilinear_vars
%     constr_back = subs(constr_back, bilinear_vars(i), prod(w{i}));
%   end
%   valuecheck(0, double(constr_back - constraint));
  
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
