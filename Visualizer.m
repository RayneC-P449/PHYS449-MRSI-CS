classdef Visualizer
   properties    
      
   end
   methods
       function obj = Visualizer()
       end

       function visualize(obj, space_v, space_h, fig_num, bg_img, vscale, hdir, vlimits)
           figure(fig_num);
           clf;
           ax_bg = axes('Units','normalized','Position',[0 0 1 1]);
           imshow(bg_img, 'Parent', ax_bg);
           axis(ax_bg, 'off');
           set(ax_bg,'DataAspectRatioMode','auto')
           voxrange = size(space_v, [1,2]);
           w = 1/voxrange(1);
           h = 1/voxrange(2);
           for r = 1:voxrange(1)
               for c = 1:voxrange(2)
                   x = (c-1)/voxrange(1);
                   y = 1-r/voxrange(2);
                   ax = axes('Units','normalized',...
                       'Position',[x y w h]);
                   plot(space_h, abs(squeeze(space_v(r,c,:)) / vscale), 'g');
                   ylim(vlimits);
                   xlim([min(space_h), max(space_h)])
                   box on;
                   ax.Color = 'None';
                   ax.XTick = [];
                   ax.YTick = [];
                   set(ax,'XDir', hdir);
               end
           end
       end
   end
end