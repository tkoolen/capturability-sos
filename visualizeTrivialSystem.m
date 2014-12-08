function visualizeTrivialSystem(B_sol, x, u_sol, f)
X = linspace(-3, 3, 100);
B_num = msubs(B_sol, x, X);
h = figure(1);
hold on;
lines = findall(gcf, 'Type', 'line');
for i = 1 : length(lines)
  set(lines(i), 'Color', 0.7 * ones(1, 3));
end
plot(X, B_num);
plot(X, zeros(size(X)), 'r');


end