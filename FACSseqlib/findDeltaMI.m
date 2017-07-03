function [deltaMI,deltapJoint]=findDeltaMI(MIlo,MIhi,pJointlo,pJointhi,name)
load('MyColormaps','mycmap')
mycmap=[[1.0 1.0 1.0];mycmap];

addpath ~/Documents/MATLAB/othercolor/othercolor/
mycmap=othercolor('RdBu5');
mycmap=mycmap(end:-1:1,:);

deltaMI=MIlo-MIhi;
deltaMI(isnan(deltaMI))=0;

lc=length(mycmap(:,1));
if max(abs(deltaMI(deltaMI<0)))<max(abs(deltaMI(deltaMI>0)))
    sc=round(max(abs(deltaMI(deltaMI<0)))/max(abs(deltaMI(deltaMI>0)))*floor(lc/2));
    mycmap=mycmap((ceil(lc/2)-sc+1):lc,:);
elseif max(abs(deltaMI(deltaMI>0)))<max(abs(deltaMI(deltaMI<0)))
    sc=floor(max(abs(deltaMI(deltaMI>0)))/max(abs(deltaMI(deltaMI<0)))*floor(lc/2));
%     max(abs(deltaMI(deltaMI>0)))
%     max(abs(deltaMI(deltaMI<0)))
%     sc
%     (floor(lc/2)+sc)
%     size(mycmap)
    mycmap=mycmap(1:(floor(lc/2)+sc),:);
end

try
    figtitle=strcat(name,': \DeltaMI');
catch
    figtitle='\DeltaMI';
end
setfig(figtitle{1});clf
% mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
imagesc(deltaMI)
c=colorbar;
ylabel(c,'\DeltaMI')
set(c,'fontsize',16)
set(c,'linewidth',1.5)
set(gca,'fontsize',20)
set(gca,'linewidth',1.5)
xlabel('position i')
ylabel('position j')
title('\DeltaMI')

deltapJoint=pJointlo-pJointhi;
addpath ~/Documents/MATLAB/othercolor/othercolor/
mycmap=othercolor('RdBu5');
mycmap=mycmap(end:-1:1,:);
lc=length(mycmap(:,1));



if max(abs(deltapJoint(deltapJoint<0)))<max(abs(deltapJoint(deltapJoint>0)))
    sc=round(max(abs(deltapJoint(deltapJoint<0)))/max(abs(deltapJoint(deltapJoint>0)))*floor(lc/2));
    sc
    mycmap=mycmap((ceil(lc/2)-sc+1):lc,:);
elseif max(abs(deltapJoint(deltapJoint>0)))<max(abs(deltapJoint(deltapJoint<0)))
    sc=round(max(abs(deltapJoint(deltapJoint>0)))/max(abs(deltapJoint(deltapJoint<0)))*floor(lc/2));
    mycmap=mycmap(1:(ceil(lc/2)+sc),:);
end


try
    figtitle=strcat(name,': \DeltaJoint probability');
catch
    figtitle='\DeltaJoint probability';
end
setfig(figtitle{1});clf
% mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
imagesc(deltapJoint)
c=colorbar;
ylabel(c,'\DeltaJoint probability')
set(c,'fontsize',20)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('\DeltaJoint probability')



end