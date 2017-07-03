% Export data from a figure into a CSV file
function exportfigdata(h,csvfile,varargin)
defaults=struct('xlabel','X','ylabel','Y','zlabel','Z');
args=processargs(defaults,varargin);

if ischar(csvfile)
  fd=fopen(csvfile,'w');
  fprintf('Saving figure data in %s\n',csvfile);
else
  fd=csvfile;
end
if isprop(h,'Title') && ~isempty(h.Title.String)
  fprintf(fd,'Title: %s\n', h.Title.String);
elseif isprop(h,'Name') && ~isempty(h.Name)
  fprintf(fd,'Name: %s\n', h.Name);
else
  fprintf(fd,'Class: %s\n', class(h));
end
if isprop(h,'XLabel') && ~isempty(h.XLabel.String)
  args.xlabel=h.XLabel.String;
end
if isprop(h,'YLabel') && ~isempty(h.YLabel.String)
  args.ylabel=h.YLabel.String;
end
if isprop(h,'ZLabel') && ~isempty(h.ZLabel.String)
  args.zlabel=h.ZLabel.String;
end
data=[];hdr={};
if isprop(h,'XData')
  fprintf('Have XData with length %d\n', length(h.XData(:)));
  data(:,end+1)=h.XData(:);
  hdr{end+1}=args.xlabel;
end
if isprop(h,'YData')
  fprintf('Have YData with length %d\n', length(h.YData(:)));
  data(:,end+1)=h.YData(:);
  hdr{end+1}=args.ylabel;
end
if isprop(h,'ZData')
  fprintf('Have ZData with length %d\n', length(h.ZData(:)));
  data(:,end+1)=h.ZData(:);
  hdr{end+1}=args.zlabel;
end
if isprop(h,'LData')
  fprintf('Have LData with length %d\n', length(h.LData(:)));
  data(:,end+1)=h.LData(:);
  hdr{end+1}='Lower';
end
if isprop(h,'UData')
  fprintf('Have UData with length %d\n', length(h.UData(:)));
  data(:,end+1)=h.UData(:);
  hdr{end+1}='Upper';
end
if ~isempty(data)
  for i=1:length(hdr)
    if i~=1
      fprintf(fd,',');
    end
    fprintf(fd,'"%s"',hdr{i});
  end
  fprintf(fd,'\n');
  for i=1:size(data,1)
    for j=1:size(data,2)
      if j~=1
        fprintf(fd,',');
      end
      fprintf(fd,'%f',data(i,j));
    end
    fprintf(fd,'\n');
  end
end
if isprop(h,'Children')
  for i=1:length(h.Children)
    exportfigdata(h.Children(i),fd,'xlabel',args.xlabel,'ylabel',args.ylabel,'zlabel',args.zlabel);
  end
end
if ischar(csvfile)
  fclose(fd);
end
