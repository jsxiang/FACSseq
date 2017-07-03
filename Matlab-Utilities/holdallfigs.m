% Used in conjuction with setfig to prevent all figure windows from being reused
function holdallfigs()
global figlist
bad=false(1,length(figlist.name));
for i=1:length(figlist.name)
  if figlist.name{i}(1)~='#'
    try
      set(figlist.fignum(i),'Name',[figlist.name{i},' (HOLD)']);
      figlist.name{i}=['#',figlist.name{i}];
    catch me
      bad(i)=true;
    end
  end
end
if any(bad)
  fprintf('Bad figures: %s\n', sprintf('%d ', find(bad)));
  figlist.name=figlist.name(~bad);
  figlist.fignum=figlist.fignum(~bad);
end

