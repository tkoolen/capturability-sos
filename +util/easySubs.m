function ret = easySubs(poly, var, input)
if isnumeric(input)
  ret = full(msubs(poly, var, input));
else
  ret = subs(poly, var, input);
end
end