function gooddata=findLengthdependence(gooddata,mus)

%% look for average length dependence of mu
m=zeros(69,69);
mc=m;
zerolengthloops=[];
for k=1:length(gooddata.seqs)
    if gooddata.loop1len(k)==0||gooddata.loop2len(k)==0
        zerolengthloops(end+1)=k;
    else
        m(gooddata.loop1len(k),gooddata.loop2len(k),end+1)=mus(k);
        mc(gooddata.loop1len(k),gooddata.loop2len(k))=mc(gooddata.loop1len(k),gooddata.loop2len(k))+1;
    end
end
msum=sum(m,3);
msum(msum==0)=NaN;
mc(mc<10)=0;
mmean=msum./mc;
setfig('loop length dependence');clf
imagesc(mmean)
load('MyColormaps','mycmap')
icmap=mycmap(end:-1:1,:);
icmap=[[1.000 1.000 1.000];icmap];
colormap(icmap)
axis([0 70 0 70])

set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 2 length')
ylabel('loop 1 length')
zlabel('\mu')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'\mu')

setfig('number per length');clf
surf(mc)
mycmap=[[1.000 1.000 1.000];mycmap];
colormap(mycmap)
% c=colorbar;
% set(c,'fontsize',14)
% set(c,'linewidth',1.5)
% ylabel(c,'count')
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 2 length')
ylabel('loop 1 length')
zlabel('no. of seqs')
axis([0 70 0 70 0 max(max(max(mc)))])

end