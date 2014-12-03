function visualize(B, Bdot, u_vertices, g_X0, g_Xu, x, f)

nmesh = 50;

[X,Y] = meshgrid(linspace(-3,3,nmesh),linspace(-3,3,nmesh));
% initial condition set
gs_X0 = full(double(msubs(g_X0,x,[X(:),Y(:)]')));
% unsafe set
gs_Xu = full(double(msubs(g_Xu,x,[X(:),Y(:)]')));
gs_B = full(double(msubs(B,x,[X(:),Y(:)]')));
gs_Bdot = cellfun(@(z) full(double(msubs(z,x,[X(:),Y(:)]'))), Bdot, 'UniformOutput', false);

[~, u_vertex_ind] = min(vertcat(gs_Bdot{:}), [], 1);
gs_u = u_vertices(u_vertex_ind);

figure(1);
hold on;
colormap(summer);
hB = surfl(X, Y, reshape(double(gs_B),nmesh,nmesh));
alpha(hB, 0.5);

[~, hB] = contour3(X, Y, reshape(gs_B, nmesh, nmesh), [0 0], 'b'); % boundary function zero level set
arrayfun(@(x) set(x, 'LineWidth', 3), hB);



r_levelset = get(hB, 'XData');
rd_levelset = get(hB, 'YData');

for j = 1 : length(r_levelset)
  x_levelset = [r_levelset{j}; rd_levelset{j}];
  path_length = [0 cumsum(sqrt(sum(diff(x_levelset, 1, 2).^2, 1)))];
  npoints = 15;
  x_levelset = interp1(path_length', x_levelset', linspace(0, path_length(end), npoints)')';
  
  ric_levelset = sum(x_levelset, 1);
  disp(['average icp location: ' num2str(mean(ric_levelset))]);
  
  f_levelset = zeros(size(x_levelset));
  for i = 1 : size(x_levelset, 2)
    x_i = x_levelset(:, i);
    Bdots = cellfun(@(z) full(double(msubs(z, x, x_i))), Bdot);
    [~, u_index] = min(Bdots);
    u_i = u_vertices(:, u_index);
    f_levelset(:, i) = f(x_i, u_i);
  end
  dt = 1e-3;
  dx_levelset = f_levelset * dt;
  dB_levelset = full(double(msubs(B, x, x_levelset + dx_levelset))) - full(double(msubs(B, x, x_levelset)));
  quiver3(x_levelset(1, :), x_levelset(2, :), zeros(1, size(x_levelset, 2)), dx_levelset(1, :), dx_levelset(2, :), dB_levelset, 'k', 'LineWidth', 3);
end

hold off;

legend(hB(1), {'B(x) = 0'});
xlabel('x1'); ylabel('x2');
title('barrier function');
view(3);
grid on;

% figure(2)
% hold on;
% % for i = 1 : numel(gs_Bdot)
% %   mesh(X, Y, reshape(double(gs_Bdot{i}),100,100)); xlabel('x1'); ylabel('x2');
% % end
% mesh(X, Y, reshape(gs_Bdot{2} - gs_Bdot{1}, nmesh, nmesh)); xlabel('x1'); ylabel('x2');
% title('barrier function derivative difference');

% figure(3);
% mesh(X, Y, reshape(gs_Bdotmin, nmesh, nmesh)); xlabel('x1'); ylabel('x2');
% title('barrier function min derivative');


% figure()
% ax = gca();
% hold on
% % contour(X,Y,reshape(gs_X0,nmesh,nmesh), [0 0],'LineWidth',3); % initial condition set
% % contour(X,Y,reshape(gs_Xu,nmesh,nmesh), [0 0],'r','LineWidth',3) % unsafe set
% % colormap(summer)
% % for i = 1 : length(u_vertices)
% %   if any(u_vertex_ind == i)
% %     c = contourf(X,Y,reshape(u_vertex_ind, nmesh, nmesh), [i i], 'b', 'LineWidth', 2); % boundary function zero level set
% %     clabel(c);
% %     u = u_vertices(i);
% %     legend_strings{end + 1} = ['u = ' num2str(u)]; %#ok<AGROW>
% %   end
% % end
% 
% % colormap(summer(size(u_vertices, 2) + 1))
% contourf(X,Y,reshape(u_vertex_ind, nmesh, nmesh), 1:size(u_vertices, 2), ''); % boundary function zero level set
% 
% % legend_strings = {};
% [c, hB] = contour(X, Y, reshape(gs_B, nmesh, nmesh), [0 0], 'b', 'LineWidth', 2); % boundary function zero level set
% 
% % ticklabels = cell(1, size(u_vertices, 2) + 1);
% % for i = 1 : size(u_vertices, 2)
% %   ticklabels{i + 1} = ['u = ' num2str(u_vertices(:, i)')];
% % end
% % cb = colorbar();
% % set(cb, 'YTick', 1 : size(u_vertices, 2));
% % set(cb, 'YTickLabel', ticklabels);
% 
% hold off;
% 
% xlabel('x');
% ylabel('dx/dt');
% 
% legend(hB, {'B(x) = 0'});
% 
% % (scaled) vector field
% % [X,Y] = meshgrid(linspace(-3,3,50),linspace(-3,3,50));
% % x1dots = reshape(double(msubs(f(1),x,[X(:),Y(:)]')),50,50);
% % x2dots = reshape(double(msubs(f(2),x,[X(:),Y(:)]')),50,50);
% % x1dots = 0.1*x1dots./(sqrt(x1dots.^2 + x2dots.^2));
% % x2dots = 0.1*x2dots./(sqrt(x1dots.^2 + x2dots.^2));
% % quiver(X,Y,x1dots,x2dots,'AutoScale','off','Color','k');


end