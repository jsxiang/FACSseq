% Save a figure as full size PDF
function pdfsavefig(fig,filename)
if nargin==1 && strcmp(fig,'all')
  global figlist;
  for i=1:length(figlist.fignum)
    try   % In case figure no longer exists
      fnum=figlist.fignum(i);
      pdfsavefig(fnum);
    catch me
    end
  end
  return;
end

if nargin<2
  name=get(fig,'Name');
  if isempty(name)
    name=sprintf('Fig%d',get(fig,'Number'));
  end
  filename=[name,'.pdf'];
  filename=strrep(filename,'/','_');
end
width=11*2;
height=8.5*2;
set(fig,'PaperPositionMode','manual','PaperUnits','inches','PaperSize',[width,height],'PaperPosition',[0 0 width height]);
set(get(gcf,'Children'),'FontSize',12);
print(fig,filename,'-dpdf','-r600');
fprintf('Saved figure %d (%s) as %s\n', get(fig,'Number'),get(fig,'Name'), filename);
set(fig,'PaperPositionMode','auto');
