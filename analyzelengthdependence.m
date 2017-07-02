addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/FACSseq/
% addpath ~/Documents/MATLAB/BREWER/
%%
% clear
g=load('YFSI_gooddata.mat');
gooddata=g.gooddata;

%% look for average length dependence of mu
m=zeros(69,69);
mc=m;
zerolengthloops=[];
for k=1:length(gooddata.seqs)
    if gooddata.loop1len(k)==0||gooddata.loop2len(k)==0
        zerolengthloops(end+1)=k;
    else
        m(gooddata.loop1len(k),gooddata.loop2len(k),end+1)=gooddata.VYBmus(k);
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
%% pick the length that has the most number of counts
i=find(mc==max(max(mc)));
l1=rem(i,69);
l2=(i-rem(i,69))/69+1;
mc(l1,l2);

%% get all the loop1 seqs that have length l1
lseqs=gooddata.loop2(gooddata.loop2len==13);

setfig('seqlogos');clf
seqlogo(lseqs)
