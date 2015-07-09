function ret = isLinearSOSConstraint(constraint, indeterminates)
[~, degrees, ~] = decomp(constraint, indeterminates);
total_degrees = sum(degrees, 2);
ret = ~any(total_degrees > 1);
end
