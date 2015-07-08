classdef SubFigureVisualizer < visualizers.ViableCaptureVisualizer
  properties (Hidden)
    hfig = [];
    cols;
    rows;
    plot_function;
    current_index = 1;
  end
  
  methods
    function obj = SubFigureVisualizer(n_subfigures, n_rows, plot_function)
      obj.cols = ceil(n_subfigures / n_rows);
      obj.rows = n_rows;
      obj.plot_function = plot_function;
    end
    
    function visualize(obj, B, x, u, f)
      if isempty(obj.hfig)
        obj.hfig = figure('Position', [100, 100, 1080, 720]);
      end
      subplot(obj.rows, obj.cols, obj.current_index);
      options.show_legend = false;
      obj.plot_function(B, x, u, f, obj.hfig, options);
      obj.current_index = obj.current_index + 1;
    end
  end
end
