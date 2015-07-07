classdef SubFigureLIPMVisualizer < ViableCaptureVisualizer
  properties (Hidden)
    hfig = [];
    cols;
    rows;
    current_index = 1;
  end
  
  methods
    function obj = SubFigureLIPMVisualizer(n_subfigures, n_rows)
      obj.cols = ceil(n_subfigures / n_rows);
      obj.rows = n_rows;
    end
    
    function visualize(obj, B, x, u, f)
      if isempty(obj.hfig)
        obj.hfig = figure('Position', [100, 100, 1080, 720]);
      end
      subplot(obj.rows, obj.cols, obj.current_index);
      options.show_legend = false;
      visualizeLIPM(B, x, u, f, obj.hfig, options);
      obj.current_index = obj.current_index + 1;
    end
  end
end
