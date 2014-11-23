function xd = lipmDynamics(x, u)
if nargin < 2
  u = 0;
end

r = x(1);
rd = x(2);
rdd = r - u;
xd = [rd; rdd];
end