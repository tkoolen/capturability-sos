function ret = functionHandleSubs(fun, varargin)
if ~iscell(fun)
  fun = {fun};
end
n = length(fun);
ret = zeros(n, 1, 'msspoly');
for i = 1 : n
  ret(i) = fun{i}(varargin{:});
end
end