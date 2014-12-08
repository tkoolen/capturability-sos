function ret = captureLimit(t_min, u_max, s_max, N)

expMinusDeltaTMin = exp(-t_min);

% N infinite
if isinf(N)
  ret = u_max + s_max * expMinusDeltaTMin / (1 - expMinusDeltaTMin);
else
  % N = 0
  ret = u_max;
  
  % TODO: use explicit formula
  % N finite
  for i = 1 : N
    ret = (s_max + ret - u_max) * expMinusDeltaTMin + u_max;
  end
end

end