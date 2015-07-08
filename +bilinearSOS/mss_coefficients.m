function [coefficients, monomials] = mss_coefficients(f, x)
[p, M, monomials] = decomp(f, x);
coefficients = recomp(p, M, speye(size(M,1)));
end