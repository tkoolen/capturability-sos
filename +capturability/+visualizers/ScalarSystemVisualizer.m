classdef ScalarSystemVisualizer < capturability.visualizers.ViableCaptureVisualizer
  properties (Access = private)
    x_min;
    x_max;
  end
  
  methods
    function obj = ScalarSystemVisualizer(x_min, x_max)
      obj.x_min = x_min;
      obj.x_max = x_max;
    end
      
    function visualize(obj, B, x, u, f)
      X = linspace(obj.x_min, obj.x_max, 100);
      B_num = msubs(B, x, X);
      figure(1);
      hold on;
      lines = findall(gcf, 'Type', 'line');
      
      for i = 1 : length(lines)
        set(lines(i), 'Color', 0.7 * ones(1, 3));
      end
%       plot(X, zeros(size(X)), 'r');
%       plot(X, B_num);
%       plot(X, msubs(f(x, u), x, X), 'g');
      plotyy(X, B_num, X, msubs(f(x, u), x, X));
      grid on;
    end
  end
end
