function visualizeScalarSystem(B_sol, x, u, f)
X = linspace(-2, 2, 100);
B_num = msubs(B_sol, x, X);
figure(1);
hold on;
lines = findall(gcf, 'Type', 'line');

for i = 1 : length(lines)
  set(lines(i), 'Color', 0.7 * ones(1, 3));
end
plot(X, B_num);
plot(X, zeros(size(X)), 'r');
plot(X, msubs(f(x, u), x, X), 'g');

end