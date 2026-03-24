classdef Visualizer
   properties    
      
   end
   methods
       function obj = Visualizer()
       end

       function visualize(~, x, Y, xdir, ylimits, img, path,toshow)
           if toshow
               f = figure('Visible', 'on');
           else
               f = figure('Visible', 'off');
           end
           ax_bg = axes('Units','normalized','Position',[0 0 1 1]);
           if numel(img) ~= 0
               imshow(img, 'Parent', ax_bg);
           end
           axis(ax_bg, 'off');
           set(ax_bg,'DataAspectRatioMode','auto')
           voxrange = size(Y, [1,2]);
           w = 1/voxrange(1);
           h = 1/voxrange(2);
           % set(gcf, 'Units', 'pixels');
           % pos = get(gcf, 'Position');
           % side = min(pos(3), pos(4));
           % set(gcf, 'Position', [pos(1), pos(2), side, side]);
           for r = 1:voxrange(1)
               for c = 1:voxrange(2)
                   rv = (c-1)/voxrange(1);
                   cv = 1-r/voxrange(2);
                   ax = axes('Units','normalized',...
                       'Position',[rv cv w h]);
                   plot(x, abs(squeeze(Y(r,c,:))), 'g');
                   ylim(ylimits);
                   xlim([min(x), max(x)])
                   box on;
                   ax.Color = 'None';
                   ax.XTick = [];
                   ax.YTick = [];
                   set(ax,'XDir', xdir);
               end
           end
           if ~toshow
               saveas(f, path);
               close(f);
           end
       end
   end
end