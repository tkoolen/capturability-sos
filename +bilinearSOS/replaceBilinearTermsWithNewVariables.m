function [prog, monomials_linear, x, X, psd_mat] = replaceBilinearTermsWithNewVariables(prog, monomials)
% replaces bilinear terms in monomials with variables from a symmetric
% matrix X, such that if we substitute X = x * x', then the original
% monomials are reobtained. Additionally, adds a PSD constraint to prog
% that 
%     [X  x;
%      x' 1]
% be positive definite, a relaxation of the condition that X = x * x'
% (dropping the constraint that this matrix should additionally be rank 1).

[vars, degrees, M] = decomp(monomials);

% check that input is indeed bilinear
total_degrees = sum(degrees, 2);
if any(total_degrees > 2) || any(any(degrees > 1)) || any(any(degrees < 0))
  error('input is not bilinear')
end

% find monomials which are bilinear and variables which appear in bilinear
% terms
bilinear_monomial_indices = find(total_degrees == 2);
degrees_bilinear = degrees(bilinear_monomial_indices, :);
x_ind = any(degrees_bilinear > 0, 1);
x = vars(x_ind);
degrees_bilinear = degrees_bilinear(:, x_ind);

% create symmetric X matrix
n_Wvec = spotprog.psdDimToNo(length(x));
[prog, Xvec] = prog.newFree(n_Wvec, 1);
X = mss_v2s(Xvec);

% create matrix that can be used to look up indices of elements of X in
% Xvec
Xvec_index_lookup = mss_v2s(1 : n_Wvec);

% find indices into Xvec corresponding to bilinear terms in monomials
Xvec_indices = zeros(size(degrees_bilinear, 1), 1);
for i = 1 : size(degrees_bilinear, 1)
  var_indices = find(degrees_bilinear(i, :));
  Xvec_indices(i) = Xvec_index_lookup(var_indices(1), var_indices(2));
end

% change degrees so that we use the bilinear variables instead of products
% of the original variables
degrees(bilinear_monomial_indices, :) = 0;
[row, col, value] = find(degrees);
row = [row; bilinear_monomial_indices];
col = [col; length(vars) + Xvec_indices];
value = [value; ones(size(Xvec_indices))];
degrees = sparse(row, col, value, size(degrees, 1), length(vars) + length(Xvec));

% recompose coefficients after getting rid of bilinear monomials and adding
% the bilinear variables instead
monomials_linear = recomp([vars; Xvec], degrees, M);

% set up the PSD constraint
psd_mat = [X x; x' 1];
prog = prog.withPSD(psd_mat);

% check
if ~isempty(x)
  X_check = x * x';
  X_checkvec = mss_s2v(X_check);
  monomials_back = subs(monomials_linear, Xvec, X_checkvec);
  valuecheck(zeros(size(monomials)), double(monomials_back - monomials));
end

end