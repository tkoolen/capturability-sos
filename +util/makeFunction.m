function fun = makeFunction(poly, var)
fun = @(input) easySubs(poly, var, input);
end