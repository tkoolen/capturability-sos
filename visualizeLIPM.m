function visualizeLIPM(B, x, u, f, videoWriter)

nmesh = 50;

[X,Y] = meshgrid(linspace(-3,3,nmesh),linspace(-3,3,nmesh));
if length(x) == 2
  gs_B = reshape(full(double(msubs(B,x,[X(:),Y(:)]'))), nmesh, nmesh);
else
  T = ones(size(X));
  gs_B = reshape(full(double(msubs(B,x,[X(:),Y(:),T(:)]'))), nmesh, nmesh);
end
% gs_Bdot = full(double(msubs(Bdot,x,[X(:),Y(:)]')));

hfig = figure('Position', [100, 100, 1080, 720]);
hold on;
colormap(summer);
hB = surfl(X, Y, reshape(double(gs_B),nmesh,nmesh));
alpha(hB, 0.5);
shading flat;
view(3);
axis tight;
grid on;
ax = gca();
zticks = get(ax, 'ZTick');
set(ax, 'ZTickMode', 'manual'); set(ax, 'ZLimMode', 'manual')
contour3(X, Y, gs_B, zticks, 'k:');
[~, hB] = contour3(X, Y, gs_B, [0 0], 'b'); % boundary function zero level set

arrayfun(@(x) set(x, 'LineWidth', 3), hB);

r_levelset = get(hB, 'XData');
rd_levelset = get(hB, 'YData');
B_levelset = get(hB, 'ZData');

if ~iscell(r_levelset)
  r_levelset = {r_levelset};
  rd_levelset = {rd_levelset};
  B_levelset = {B_levelset};
end;

for j = 1 : length(r_levelset)
  if length(x) == 2
    x_levelset_j = [r_levelset{j}; rd_levelset{j}];
  else
    x_levelset_j = [r_levelset{j}; rd_levelset{j}; zeros(1, size(r_levelset{j}, 2))];
  end

  B_levelset_j = B_levelset{j};
  
  if ~isempty(B_levelset_j)
    path_length = [0 cumsum(sqrt(sum(diff(x_levelset_j, 1, 2).^2, 1)))];
    points_per_length_unit = 2;
    npoints = ceil(path_length(end) * points_per_length_unit);
    x_levelset_j = interp1(path_length', x_levelset_j', linspace(0, path_length(end), npoints)')';
    B_levelset_j = interp1(path_length', B_levelset_j', linspace(0, path_length(end), npoints)')';
    
    ric_levelset = sum(x_levelset_j, 1);
    disp(['max icp distance: ' num2str(max(abs(ric_levelset)))]);
    
%     Bdot_levelset = full(double(msubs(Bdot, x, x_levelset_j)));
    u_levelset_j = full(double(msubs(u, x, x_levelset_j)));
    
    f_levelset = zeros(size(x_levelset_j));
    for i = 1 : size(x_levelset_j, 2)
      x_i = x_levelset_j(:, i);
      u_i = u_levelset_j(:, i);
      f_levelset(:, i) = f(x_i, u_i);
    end
    dt = 1e-3;
    dx_levelset = f_levelset * dt;
    dB_levelset = full(double(msubs(B, x, x_levelset_j + dx_levelset))) - full(double(msubs(B, x, x_levelset_j)));
    
    quiver3(x_levelset_j(1, :), x_levelset_j(2, :), B_levelset_j, dx_levelset(1, :), dx_levelset(2, :), dB_levelset, 'r', 'LineWidth', 2);
  end
end

hold off;

set(gca,'FontSize', 15)
legend(hB(1), {'B(x) = 0'});
axis([-3 3 -3 3 -5 35]);
xlabel('r'); ylabel('dr/dt');

if nargin > 4
  frame = getframe(hfig);
  writeVideo(videoWriter,frame);
end

end