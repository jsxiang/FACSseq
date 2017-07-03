% Save a figure for use in a publication
% Sets width, height in inches
function pubfigure(filename,fnum,width,height,varargin)
  defaults=struct('tickfontsize',10,'axisfontsize',12,'titlefontsize',14,'format','epsc2','markersize',[],'scale',1,'savedata',false);
  args=processargs(defaults,varargin);
  if nargin<2
    fnum=gcf;
  end
  if nargin<3
    width=6;
  end
  if nargin<4
    height=width*2/3;
  end
  width=width*args.scale;
  height=height*args.scale;
  figure(fnum);

  %  set(gca,'Box','off');
  % left, bottom, right, top margins
  ps=get(gca,'Position'); % Need to handle stuff at right, such as colorbars
  margins=args.axisfontsize*[3.5 3 1 1]/72;
  if ~isempty(get(get(gca,'Title'),'String'))
    margins(4)=args.titlefontsize*2.5/72;
  end
  ps(1)=margins(1)/width;
  ps(2)=margins(2)/height;
  ps(4)=1-(margins(2)+margins(4))/height;
  % set(gca,'Position',ps);
  set(gcf,'PaperUnits','inches');
  set(gcf,'PaperPosition',[0 0 width height]);
  set(gcf,'PaperSize',[width height]);
  set(gcf,'Units','inches');
  ps=get(gcf,'Position');
  set(gcf,'Position',[ps(1) ps(2) width height]);
  axes=get(gcf,'children');
  for j=1:length(axes)
    if ~strcmp(class(axes(j)),'matlab.graphics.axis.Axes')
      continue;
    end
    set(axes(j),'FontSize',args.tickfontsize);
    set(get(axes(j),'XLabel'),'FontSize',args.axisfontsize,'FontWeight','normal');
    set(get(axes(j),'YLabel'),'FontSize',args.axisfontsize,'FontWeight','normal');
    set(get(axes(j),'Title'),'FontSize',args.titlefontsize,'FontWeight','bold');
    c=get(axes(j),'Children');
    if ~isempty(args.markersize)
      for i=1:length(c)
        try
          if get(c(i),'Marker')=='x'
            set(c(i),'MarkerSize',8);
          else
            set(c(i),'MarkerSize',args.markersize);
          end
        catch me
          %        fprintf('Ignoring exception\n');
        end
      end
    end
  end
  print(gcf,sprintf('-d%s',args.format),filename);
  fprintf('Saved figure to %s with format %s\n', filename,args.format);
  if args.savedata
    % Also export data
    exportfigdata(gcf,[filename,'.csv']);
  end
end
