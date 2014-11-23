function frlibProblem()
% generated problem using:

% prog = spotsosprog;
% [prog, x] = prog.newIndeterminate('x', 1);
% [prog, V] = prog.newSOSPoly(monomials(x, 0 : 4));
% solver = @spot_mosek;
% options = spot_sdp_default_options();
% options.verbose = 1;
% prog.minimize(0,solver,options);

load('frlib_problem.mat');
use_frlib = true;
if use_frlib
  frPrg = frlibPrg(A,b,c,K);
%   opts.useQR = 1;
%   frPrg = frPrg.ReducePrimal('d', opts);
  A = frPrg.A;
  b = frPrg.b;
  c = frPrg.c';
  K = frPrg.K;
end

options = spot_sdp_default_options();
options.verbose = 1;
spot_mosek(A,b,c,K,options);
end