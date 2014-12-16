function testViableCaptureContinuous(verify_manual)
if nargin < 1
  verify_manual = true;
end
% path setup
addpath(fullfile('util'));
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
g_Xg = @(x) -(x + 1.9)^2 / 0.1^2  + 1;
g_Xf = @(x) x^2 / 2^2 - 1;

options = struct;
options.plotfun = @visualizeScalarSystem;
if verify_manual
%   options.B_manual = @(x) x^2 / 1^2 - 1; %feasible with X_g = -1
%   options.B_manual = @(x) x^2 / 2^2 - 1; %feasible with X_g = 1
  options.B_manual = @(x) (x + 0.5)^2 / 1.49^2 - 1; % feasible with X_g positive only in one corner
end
[B_fun, u_fun] = viableCaptureContinuous(f, nstates, u_min, u_max, g_Xg, g_Xf, options);
end