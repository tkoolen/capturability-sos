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
valuecheck(zeros(size(cconstr)), double(subs(cconstr_W, mss_s2v(W), mss_s2v(ww)) - cconstr));
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