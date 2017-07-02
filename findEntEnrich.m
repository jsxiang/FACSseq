function [deltaEnt,deltaPAll]=findEntEnrich(ent25,ent75,pAll25,pAll75,name)
% find difference between entropy
H0=pAll75(1,:).*log2(pAll75(1,:))+pAll75(2,:).*log2(pAll75(2,:))+pAll75(3,:).*log2(pAll75(3,:))+pAll75(3,:).*log2(pAll75(3,:));

H1=pAll25(1,:).*log2(pAll25(1,:))+pAll25(2,:).*log2(pAll25(2,:))+pAll25(3,:).*log2(pAll25(3,:))+pAll25(3,:).*log2(pAll25(3,:));


deltaEnt=H1-H0;
deltaPAll=pAll25-pAll75;
addpath ~/Documents/MATLAB/cbrewer/cbrewer/cbrewer
mycmap=cbrewer('qual','Paired',7);
colormap(mycmap)

try
    figtitle=strcat(name,': Entropy enrichment');
catch
    figtitle='Entropy enrichment';
end
setfig(figtitle{1});clf


% bar([deltaPAll.*repmat(deltaEnt,4,1)]','Stacked','FaceAlpha',0.75)
% legend('A','U','C','G','Location','Best')
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% ylabel('entropy')
% xlabel('position')

colormap(mycmap)

deltaEnt(deltaEnt<0)=0;
entforbar=[deltaPAll.*repmat(deltaEnt,4,1)]';
% minent=min(deltaEnt);
% entforbar=[deltaPAll.*repmat(deltaEnt-minent,4,1)]';
Xneg=entforbar;
Xpos=entforbar;
Xneg(entforbar>0) = 0;
Xpos(entforbar<0) = 0;
hold on
bar(Xneg,'stack')
bar(Xpos,'stack')

hold off

% bar(entforbar,'Stacked')
% colormap('jet')
legend('A','U','C','G','Location','Best')
set(gca,'fontsize',20)
set(gca,'linewidth',1.5)
ylabel('\Deltaentropy')
xlabel('position i')
end