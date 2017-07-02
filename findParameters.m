function findParameters(mus,seqs)
% function determines the different percentile mus and displays them for
% visual inspection 
% input requires mu values previously determined, and aligned sequences

l1l2seqs=seqs;
totallooplength=length(l1l2seqs(1,:));
setfig('check distribution');clf
DNA={'A','T','C','G'};
for p=1:totallooplength
    for i=1:length(DNA)
        subplot(totallooplength,4,4*(p-1)+i)
        ind=l1l2seqs(:,p)==i;
        m=mus(ind);
        hist(m,12);
        title(DNA{i})
        xlabel('\mu')
        ylabel('counts')
        xlim([0 12])
%         set(gca,'fontsize',14)
%         set(gca,'linewidth',1.5)
        p5(p,i)=prctile(m,5);
        p10(p,i)=prctile(m,10);
        p25(p,i)=prctile(m,25);
        md(p,i)=median(m);
    end
end

setfig('find parameters');clf
subplot(1,4,1)
imagesc(p5)
c=colorbar;
set(c,'fontsize',12)
set(c,'linewidth',1.5)
ylabel(c,'\mu')

load('MyColormaps','mycmap')
icmap=mycmap;
icmap=[[1.000 1.000 1.000];icmap];
colormap(icmap)
L = get(gca,'XLim');
xlabelnames=DNA;
NumTicks = length(xlabelnames)+5;
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{'','A','','T','','C','','G',''})
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('5th percentile')

subplot(1,4,2)
imagesc(p10)
c=colorbar;
set(c,'fontsize',12)
set(c,'linewidth',1.5)
ylabel(c,'\mu')

colormap(icmap)

L = get(gca,'XLim');
xlabelnames=DNA;
NumTicks = length(xlabelnames)+5;
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{'','A','','T','','C','','G',''})
title('10th percentile')

subplot(1,4,3)
imagesc(p25)
c=colorbar;
set(c,'fontsize',12)
set(c,'linewidth',1.5)
ylabel(c,'\mu')

colormap(icmap)
L = get(gca,'XLim');
xlabelnames=DNA;
NumTicks = length(xlabelnames)+5;
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{'','A','','T','','C','','G',''})
title('25th percentile')

subplot(1,4,4)
imagesc(md)
c=colorbar;
set(c,'fontsize',12)
set(c,'linewidth',1.5)
ylabel(c,'\mu')

colormap(mycmap)

L = get(gca,'XLim');
xlabelnames=DNA;
NumTicks = length(xlabelnames)+5;
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{'','A','','T','','C','','G',''})
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('50th percentile')




end