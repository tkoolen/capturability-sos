classdef SeparateFrameLIPMVisualizer < LIPMVisualizer
  methods
    function visualize(obj, B, x, u, f)
      hfig = figure('Position', [100, 100, 1080, 720]);
      options.show_legend = true;
      visualizeLIPM(B, x, u, f, hfig, options);
    end
  end
end