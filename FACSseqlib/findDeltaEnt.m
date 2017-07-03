function [deltaEnt,deltaPAll]=findDeltaEnt(ent25,ent75,pAll25,pAll75,name)
% find difference between entropy
deltaEnt=ent25-ent75;
deltaPAll=pAll25-pAll75;
addpath ~/Documents/MATLAB/cbrewer/cbrewer/cbrewer
mycmap=cbrewer('qual','Set2',40);
colormap(mycmap)

try
    figtitle=strcat(name,': \Deltaentropy');
catch
    figtitle='\Deltaentropy';
end
setfig(figtitle{1});clf

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