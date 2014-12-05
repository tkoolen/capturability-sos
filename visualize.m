function visualize(B, Bdots, u_vertices, x, u, f)

nmesh = 50;

[X,Y] = meshgrid(linspace(-3,3,nmesh),linspace(-3,3,nmesh));
% initial condition set
% gs_X0 = full(double(msubs(g_X0,x,[X(:),Y(:)]')));
% unsafe set
% gs_Xu = full(double(msubs(g_Xu,x,[X(:),Y(:)]')));
gs_B = full(double(msubs(B,x,[X(:),Y(:)]')));
gs_Bdot = cellfun(@(z) full(double(msubs(z,x,[X(:),Y(:)]'))), Bdots, 'UniformOutput', false);

[~, u_vertex_ind] = min(vertcat(gs_Bdot{:}), [], 1);
% gs_u = u_vertices(u_vertex_ind);

figure();
hold on;
colormap(summer);
hB = surfl(X, Y, reshape(double(gs_B),nmesh,nmesh));
alpha(hB, 0.5);

[~, hB] = contour3(X, Y, reshape(gs_B, nmesh, nmesh), [0 0], 'b'); % boundary function zero level set
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
  x_levelset_j = [r_levelset{j}; rd_levelset{j}];
  B_levelset_j = B_levelset{j};
  
  if ~isempty(B_levelset_j)
    path_length = [0 cumsum(sqrt(sum(diff(x_levelset_j, 1, 2).^2, 1)))];
    points_per_length_unit = 2;
    npoints = ceil(path_length(end) * points_per_length_unit);
    x_levelset_j = interp1(path_length', x_levelset_j', linspace(0, path_length(end), npoints)')';
    B_levelset_j = interp1(path_length', B_levelset_j', linspace(0, path_length(end), npoints)')';
    
    ric_levelset = sum(x_levelset_j, 1);
    disp(['max icp distance: ' num2str(max(abs(ric_levelset)))]);
    
    Bdots_levelset = cellfun(@(z) full(double(msubs(z, x, x_levelset_j))), Bdots, 'UniformOutput', false);
    Bdots_levelset = vertcat(Bdots_levelset{:});
    [~, u_indices_levelset] = min(Bdots_levelset, [], 1);
    u_levelset_j = u_vertices(:, u_indices_levelset);
    
    f_levelset = zeros(size(x_levelset_j));
    for i = 1 : size(x_levelset_j, 2)
      x_i = x_levelset_j(:, i);
      u_i = u_levelset_j(:, i);
      f_levelset(:, i) = f(x_i, u_i);
    end
    dt = 1e-3;
    dx_levelset = f_levelset * dt;
    dB_levelset = full(double(msubs(B, x, x_levelset_j + dx_levelset))) - full(double(msubs(B, x, x_levelset_j)));
    
    hU = zeros(size(u_vertices, 2), 1);
    u_legend_strings = cell(size(u_vertices, 2), 1);
    for i = 1 : size(u_vertices, 2)
      u_i_indices = u_indices_levelset == i;
      hU(i) = quiver3(x_levelset_j(1, u_i_indices), x_levelset_j(2, u_i_indices), B_levelset_j(u_i_indices), dx_levelset(1, u_i_indices), dx_levelset(2, u_i_indices), dB_levelset(u_i_indices), 'LineWidth', 2);
      u_legend_strings{i} = ['u = ' mat2str(u_vertices(:, i)', 2)];
    end
  end
end

hold off;

legend([hB(1); hU], {'B(x) = 0', u_legend_strings{:}});
xlabel('r'); ylabel('dr/dt');
view(3);
grid on;


end