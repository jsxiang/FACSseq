% Add logarithmic ticks to an axis
function logticks(xlog,ylog)
if nargin<1
  xlog=true;
end
if nargin<2
  ylog=xlog;
end
c=axis;
for ax=1:2
  if ax==1 && ~xlog
    continue;
  end
  if ax==2 && ~ylog
    continue;
  end
  ticks=[];
  ticklabels={};
  if c(ax*2-1)<=0
    if c(ax)<=0
      error('Can''t add log ticks since axes bounds are <= 0');
    end
    c(ax*2-1)=c(ax)/1000;
  end
  ndecades=log10(c(ax*2))-log10(c(ax*2-1));
  for i=floor(log10(c(ax*2-1))):ceil(log10(c(ax*2)))
    if i<0
      fmt=sprintf('%%.%df',-i);
    elseif i>=4
      fmt='%.1g';
    else
      fmt='%.0f';
    end
    for j=1:9
      tval=j*10^i;
      if tval<c(ax*2-1) || tval>c(ax*2)
        continue;
      end
      ticks(end+1)=tval;
      if j==1 || (j==2 && ndecades<5) || (j==5 && ndecades <3)
        if abs(tval)>=1 && abs(tval)<10000
          ticklabels{end+1}=sprintf('%d',tval);
        elseif abs(tval)>=0.1 && abs(tval)<1
          ticklabels{end+1}=sprintf('%.1f',tval);
        elseif tval/10^i == 1
          ticklabels{end+1}=sprintf('10^{%d}',i);
        else
          ticklabels{end+1}=sprintf('%dx10^{%d}',tval/10^i,i);
        end
      else
        ticklabels{end+1}='';
      end
    end
  end
  if ax==1
    set(gca,'XTick',ticks);
    set(gca,'XTickLabel',ticklabels);
  else
    set(gca,'YTick',ticks);
    set(gca,'YTickLabel',ticklabels);
  end
end
