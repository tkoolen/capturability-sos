function lipmContinuousBarrier(N, verify_manual_barrier_function)
% path setup
import examples.lipm.*
import util.*
import capturability.*;
import capturability.visualizers.*;

checkDependency('spotless');
oldpath = cd(fullfile('..', 'frlib'));
setup();
cd(oldpath);

% parameters
if nargin < 1
  N = 0;
end
if nargin < 2
  verify_manual_barrier_function = true;
end

s_max = 1;

% dynamics
f = @(x, u) [lipmDynamics(x, u); 1];

% input limits
u_min = -1;
u_max = 1;

% Guard
t_min = 1;
g_Xguard = @(x) x(3) - t_min;

if N > 0
  dN_minus_one = captureLimit(t_min, u_max, s_max, N - 1);
  BN_minus_one_prime = @(x) (x(1) + x(2))^2 / (dN_minus_one + s_max)^2 - 1;
  g_Xtarget = {@(x) -BN_minus_one_prime(x); g_Xguard};
else
  g_Xtarget = {};
end

nstates = 3;

% failed states
x_ic_dist = 3;
g_Xf = @(x) (x(1) + x(2))^2 - x_ic_dist^2;

% initial states
x_star = [0; 0];
g_Xstar = @(x) -(x(1:2) - x_star)' * (x(1:2) - x_star);

% visualization
if verify_manual_barrier_function
  options.visualizer = SeparateFrameVisualizer(@visualizeLIPM);
else
  % options.visualizer = SeparateFrameVisualizer(@visualizeLIPM);
  % options.visualizer = VideoVisualizer(@visualizeLIPM);
  options.visualizer = SubFigureVisualizer(15, 3, @visualizeLIPM);
end


margin = 0;
if verify_manual_barrier_function
  if N == 0
    dN = captureLimit(t_min, u_max, s_max, N);
    options.B_manual = @(x) (x(1) + x(2))^2 / (dN - margin)^2 - 1;
%     options.B_manual = @(x) x(1:2)' * x(1:2) / 5e-1^2 - 1;
  else
    dN_minus_one = captureLimit(t_min, u_max, s_max, N - 1);
    options.B_manual = @(x) ((x(1) + x(2))^2 / (dN_minus_one + s_max - margin)^2 - 1);
  end
end

[BN_fun, u_fun] = viableCaptureContinuous(f, nstates, u_min, u_max, g_Xtarget, g_Xf, g_Xstar, options);

end