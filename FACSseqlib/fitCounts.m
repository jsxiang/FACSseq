function [m,s,sim]=fitCounts(counts,edges,pltfig)
cr=round(counts);
cOr=counts./cr;
sim=zeros(1,sum(cr));
er=cOr.*edges;
sim(1:cr(1))=ones(1,cr(1))*er(1);
for i=2:length(cr)
    startind=sum(cr(1:i-1))+1;
    endind=sum(cr(1:i-1))+cr(i);
    sim(startind:endind)=ones(1,cr(i))*er(i);
end


[m,s]=normfit(sim);


if pltfig
    plot(er,cr,'o','MarkerSize',10,'LineWidth',1.5)
    hold on
    xsim=linspace(1,length(cr));
    ysim=normpdf(xsim,m,s)*sum(cr);
    plot(xsim,ysim,'LineWidth',1.5);
    pt=sprintf('\\mu = %s', num2str(m));
    text(m*1.05,normpdf(m,m,s)*sum(cr)*1.01,pt)

    xlabel('Bin #')
    ylabel('Count')
    legend('FACS-seq data','Norm fit','Location','NorthOutside')
    set(gca,'Linewidth',1.5)
    set(gca,'FontSize',13)
    hold off
end
end

