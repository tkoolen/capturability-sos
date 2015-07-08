function testViableCaptureContinuous(verify_manual)
if nargin < 1
  verify_manual = true;
end

% path setup
import capturability.*;
import capturability.visualizers.*;
checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

nstates = 1;
f = @(x, u) x * (x - 1) * (x + 1) + u;
u_min = 0;
u_max = 0;
% g_Xg = @(x) 1;
% g_Xg = @(x) -1;
% g_Xg = @(x) -(x + 1.9)^2 / 0.1^2  + 1; % -2 to -1.8
g_Xg = @(x) -1;
g_Xf = @(x) x^2 / 3^2 - 1; % outside [-3, 3]

options = struct;
options.visualizer = ScalarSystemVisualizer(-3.5, 3.5);
if verify_manual
%   options.B_manual = @(x) x^2 / 1^2 - 1; %feasible with X_g = -1
%   options.B_manual = @(x) x^2 / 2^2 - 1; %feasible with X_g = 1
  options.B_manual = @(x) (x + 0.5)^2 / 1.49^2 - 1; % feasible with X_g positive only in one corner
end
g_Xstar = @(x) -x^2;
[B_fun, u_fun] = viableCaptureContinuous(f, nstates, u_min, u_max, g_Xg, g_Xf, g_Xstar, options);
end