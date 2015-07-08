function [] = lipmExactPlot()
import examples.lipm.*

nmesh = 50;
rmax = 3;
rdmax = 3;
r = linspace(-rmax, rmax, nmesh);
rd = linspace(-rdmax, rdmax, nmesh);
[R, RD] = meshgrid(r, rd);

Nmax = 3;
Ns = [(0 : Nmax)'; inf];
d = zeros(size(Ns));
labels = cell(1, size(Ns, 1));
t_min = 1;
u_max = 1;
s_max = 1;
for i = 1 : length(Ns)
  N = Ns(i);
  d(i) = captureLimit(t_min, u_max, s_max, N);
  labels{i} = ['{\itN} = ' num2str(N)];
end

RIC = abs(R + RD);
figure();
colormap(summer)
[~, hContourGroup] = contourf(R, RD, -RIC, -d);
hContours = get(hContourGroup, 'children');
legend(hContours, labels);
axis equal
axis tight;
xlabel('r'); ylabel('dr/dt');
grid on;
exportFigure('../../report/figures/lipmExact', 250, 250);

end

