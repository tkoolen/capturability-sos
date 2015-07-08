function basis = monomialSubs(name, vars, degrees)
x = msspoly(name, length(vars));
basis = monomials(x, degrees);
basis = subs(basis, x, vars);
end