classdef SeparateFrameVisualizer < capturability.visualizers.ViableCaptureVisualizer
  properties (Hidden)
    plot_function;
  end
  
  methods
    function obj = SeparateFrameVisualizer(plot_function)
      obj.plot_function = plot_function;
    end
    
    function visualize(obj, B, x, u, f)
      hfig = figure('Position', [100, 100, 1080, 720]);
      options.show_legend = true;
      obj.plot_function(B, x, u, f, hfig, options);
    end
  end
end