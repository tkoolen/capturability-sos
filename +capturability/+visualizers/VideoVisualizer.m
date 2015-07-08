classdef VideoVisualizer < capturability.visualizers.ViableCaptureVisualizer
  properties (Hidden)
    hfig = [];
    videoWriter;
    plot_function;
  end
  
  methods
    function obj = VideoVisualizer(plot_function)
      writerObj = VideoWriter('barrierGrow.mp4', 'MPEG-4');
      writerObj.FrameRate = 5;
      open(writerObj);
      obj.videoWriter = writerObj;
      obj.plot_function = plot_function;
    end
    
    function delete(obj)
      close(obj.videoWriter);
    end
    
    function visualize(obj, B, x, u, f)
      if isempty(obj.hfig)
        obj.hfig = figure('Position', [100, 100, 1080, 720]);
      end
      clf(obj.hfig);
      options.show_legend = false;
      obj.plot_function(B, x, u, f, obj.hfig, options);
      drawnow;
      frame = getframe(obj.hfig);
      writeVideo(obj.videoWriter,frame);
    end
  end
end