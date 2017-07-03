% holdfig
% Used in conjuction with savefig to prevent a figure window from being reused
function holdfig(fnum)
global figlist
if ~isnumeric(fnum)
  fnum=get(fnum,'Number');
end
sel=[figlist.fignum]==fnum;
if isempty(sel)
  figlist.name{end+1}=sprintf('HOLD%d',fnum);
  figlist.fignum{end+1}=fnum;
else
  fprintf('Held figure %d: %s\n', fnum, figlist.name{sel});
  figlist.name{sel}=[figlist.name{sel},'-',sprintf('HOLD%d',fnum)];
end

