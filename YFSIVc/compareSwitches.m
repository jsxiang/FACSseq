addpath ../MATLAB/Matlab-Utilities/
addpath ../MATLAB/FACS
addpath ../MATLAB/FACSseq/
clear

y=load('yfsIVc_theo.mat');
ydata=y.theo;
y=load('yfsIVc_FA.mat');
ydata(2)=y.FA;
y=load('yfsIVc_rbz.mat');
ydata(3)=y.rbz;

m=load('mfsIV_theo.mat');
mdata=m.theo;
mdata(1).VYBmus=[mean(mdata(1).VYBmus(:,[1 3]),2) mdata(1).VYBmus(:,2)];

m=load('mfsIV_FA.mat');
mdata(2).mus=m.FA.mus;
mdata(2).sigma=m.FA.sigma;
mdata(2).seqs=m.FA.seqs;
mdata(2).combinedmus=m.FA.combinedmus;
mdata(2).combinedsigma=m.FA.combinedsigma;
mdata(2).VYBmus=[mean(m.FA.VYBmus(:,[1 3]),2) m.FA.VYBmus(:,4)];
mdata(2).origbincounts=m.FA.origbincounts;
mdata(2).scaledbincounts=m.FA.scaledbincounts;
mdata(2).rep1ind=m.FA.rep1ind;
mdata(2).rep2ind=m.FA.rep2ind;
mdata(2).rep3ind=m.FA.rep3ind;
mdata(2).rep4ind=m.FA.rep4ind;


m=load('mfsIV_rbz.mat');
mdata(3).mus=m.rbz.mus;
mdata(3).sigma=m.rbz.sigma;
mdata(3).seqs=m.rbz.seqs;
mdata(3).combinedmus=m.rbz.combinedmus;
mdata(3).combinedsigma=m.rbz.combinedsigma;
mdata(3).VYBmus=[mean(m.rbz.VYBmus(:,[1 3]),2) m.rbz.VYBmus(:,3)];
mdata(3).origbincounts=m.rbz.origbincounts;
mdata(3).scaledbincounts=m.rbz.scaledbincounts;
mdata(3).rep1ind=m.rbz.rep1ind;
mdata(3).rep2ind=m.rbz.rep2ind;
mdata(3).rep3ind=m.rbz.rep3ind;
mdata(3).rep4ind=m.rbz.rep4ind;


%% group by
common=struct;
common.seqs={};
common.yVYBmus1=[];
common.yVYBmus2=[];
common.mVYBmus1=[];
common.mVYBmus2=[];
common.origbincounts=zeros(1,length(ydata(1).origbincounts(1,:))+length(mdata(1).origbincounts(1,:))); % horizontal concatenation
common.mus=zeros(1,2);
common.sigma=zeros(1,2);

common(2)=common(1);
common(3)=common(1);

for k=1:length(mdata)
for i=1:length(mdata(k).seqs)
    iy=containers.Map(ydata(k).seqs,1:length(ydata(k).seqs));
    im=containers.Map(mdata(k).seqs,1:length(mdata(k).seqs));
    try
    common(k).yVYBmus1(end+1)=ydata(k).VYBmus(iy(mdata(k).seqs{i}),1);
    common(k).yVYBmus2(end+1)=ydata(k).VYBmus(iy(mdata(k).seqs{i}),2);
    common(k).mVYBmus1(end+1)=mdata(k).VYBmus(i,1);
    common(k).mVYBmus2(end+1)=mdata(k).VYBmus(i,2);
    common(k).seqs{end+1}=mdata(k).seqs{i};
    common(k).mus=[common(k).mus;[mdata(k).combinedmus(i) ydata(k).mus(iy(mdata(k).seqs{i}),1)]];
    common(k).sigma=[common(k).sigma;[mdata(k).combinedsigma(i) ydata(k).sigma(iy(mdata(k).seqs{i}),1)]];
    common(k).rep1ind=mdata(k).rep1ind;
    common(k).rep2ind=ydata(k).rep1ind+length(mdata(k).origbincounts(1,:));
    common(k).origbincounts=[common(k).origbincounts;[round(mdata(k).scaledbincounts(i,:),0) ydata(k).origbincounts(iy(mdata(k).seqs{i}),:)]];
    
    catch
    end
end
common(k).mus=common(k).mus(2:end,:);
common(k).sigma=common(k).sigma(2:end,:);
common(k).origbincounts=common(k).origbincounts(2:end,:);
end
%%
save('common.mat','common');
%%
lib={'Theophylline','Folinic acid','Ribozymes'};
for k=1:length(lib)
setfig(lib{k});clf
if k==3
    plot(common(k).mVYBmus1',mean([common(k).yVYBmus1' common(k).yVYBmus2'],2),'.','MarkerSize',20)
else
plot(common(k).mVYBmus1,common(k).yVYBmus1,'.','MarkerSize',20)
hold on
plot(common(k).mVYBmus2,common(k).yVYBmus2,'.','MarkerSize',20)
end

x=linspace(min(common(k).mVYBmus1),max(common(k).mVYBmus1));
hold on
mdl = fitlm(common(k).mVYBmus1,common(k).yVYBmus1);
rsq=mdl.Rsquared.Adjusted;
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
hold off
xlabel('mammalian cells')
ylabel('yeast')
set(gca,'linewidth',2)
set(gca,'fontsize',28)
t=sprintf('R^2 = %0.2f\ny = %0.2fx + %0.2f',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
if k==3
    text(0.1,0,t,'fontsize',28)
else
    text(0.4,0,t,'fontsize',28)
end
title(lib{k})



end

setfig('residuals')
plotResiduals(mdl,'probability')

setfig('diagnostic')
plotDiagnostics(mdl,'contour')

%% check controls
% spiked in controls
sTRSV='GCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGC';
g814='GCTGTCACCGGATGTTTCCGGTCTGATGAGTCCATAAGGACGAAACAGC';
g833='GCTGTCACCGGATACTTCCGGTCTGATGAGTCCCAGAGGACGAAACAGC';
g862='GCTGTCACCGGATGCATCCGGTCTGATGAGTCCCGCGGGACGAAACAGC';
L2b9a1='GCTGTCACCGGAATCAAGGTCCGGTCTGATGAGTCCGTTGTCCAATACCAGCATCGTCTTGATGCCCTTGGCAGTGGATGGGGACGGAGGACGAAACAGC';
Rs2='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAAAGGACGAAACAGC';
Theo1041='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAGAGGACGAAACAGC';
sTRSVctl='GCTGTCACCGGATGTGCTTTCCGGTACGTGAGGTCCGTGAGGACGAAACAGC';

theolibcontrols={sTRSV,g814,g833,g862,L2b9a1,Theo1041,Rs2};
ctrlnames={'sTRSV','g814','g833','g862','L2b9a1','Theo1041','Rs2'};
controlsMFSIV=zeros(length(theolibcontrols),4);

hold on
for i=1:length(common) 
        setfig(lib{i})
    for k=1:length(theolibcontrols)

    s=regexp(common(i).seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    
    ctrlmus=[common(i).mVYBmus1(ctrl)' common(i).mVYBmus2(ctrl)' common(i).yVYBmus1(ctrl)' common(i).yVYBmus2(ctrl)'];
    try
        controlsMFSIV(k,:)=ctrlmus;
    catch
    end
    hold on
    if i==3
        plot(ctrlmus(:,1),mean(ctrlmus(:,[3 4]),2),'o','linewidth',2,'markersize',16,'color',[0.8 0.4 0.1])
    else
    plot(ctrlmus(:,1),ctrlmus(:,3),'o','linewidth',2,'markersize',16,'color',[0.8 0.4 0.1])
    plot(ctrlmus(:,2),ctrlmus(:,4),'o','linewidth',2,'markersize',16,'color',[0.8 0.4 0.1])
    legend('no ligand','with ligand','location','best')
    end
    end
end

controlsMFSIV(controlsMFSIV==0)=nan;
%%
setfig('controls yeast vs mammalian');clf

% plot(controlsMFSIV(:,1),[mean(controlsMFSIV(1:3,3:4),2); controlsMFSIV(4:7,3)],'.','MarkerSize',15)
plot(controlsMFSIV(:,1),controlsMFSIV(:,3),'.','MarkerSize',16)
mdl = fitlm(controlsMFSIV(:,1),controlsMFSIV(:,3));
rsq=mdl.Rsquared.Adjusted;
xlabel('mammalian cells')
ylabel('yeast')
title('Spiked in controls')
set(gca,'linewidth',2)
set(gca,'fontsize',28)

x=linspace(min(controlsMFSIV(:,1)),max(controlsMFSIV(:,1)));
hold on
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
t=sprintf('R^2 = %0.2f\n y = %0.2fx + %0.2f',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
text(-0.25,0.1,t,'fontsize',28)
% legend('Library','1:1',ctrlnames{:},'Location','Best')
% 
% 

%% all datasets
ymu1=[common(1).yVYBmus1 common(2).yVYBmus1 common(3).yVYBmus1];
ymu2=[common(1).yVYBmus2 common(2).yVYBmus2 common(3).yVYBmus2];
mmu1=[common(1).mVYBmus1 common(2).mVYBmus1 common(3).mVYBmus1];
mmu2=[common(1).mVYBmus2 common(2).mVYBmus2 common(3).mVYBmus2];

setfig('all compare');clf
plot(mmu1,ymu1,'.','Markersize',20)
mdl = fitlm(mmu1,ymu1);
rsq=mdl.Rsquared.Adjusted;
x=linspace(min(mmu1),max(mmu1));
hold on
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
t=sprintf('R^2 = %0.2f\ny = %0.2fx + %0.2f',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
text(-0.4,2.5,t,'fontsize',28)
for i=1:length(controlsMFSIV(:,1))
%     if i<5
%     plot(controlsMFSIV(i,1),mean(controlsMFSIV(i,[3 4])),'o','color',[0.8 0.4 0.1],'linewidth',2)
%     else
    plot(controlsMFSIV(i,1),mean(controlsMFSIV(i,3)),'o','color',[0.8 0.4 0.1],'linewidth',2)
%     end
end

title('All common sequences')
xlabel('mammalian cells')
ylabel('yeast') 
set(gca,'linewidth',2)
set(gca,'fontsize',28)

%% combined then estimated variables
mus=[common(1).mus;common(2).mus;common(3).mus];
sigmas=[common(1).sigma;common(2).sigma;common(3).sigma];

setfig('morestats');clf
hold on
% herrorbar(mus(:,1),mus(:,2),sigmas(:,1),'o')
% errorbar(mus(:,1),mus(:,2),sigmas(:,2),'o')
h=herrorbar(common(1).mus(:,1),common(1).mus(:,2),common(1).sigma(:,1)./sqrt(sum(common(1).origbincounts(:,common(1).rep1ind),2)),common(1).sigma(:,1)./sqrt(sum(common(1).origbincounts(:,common(1).rep1ind),2)),'.',[0.6 0.6 0.6]);
e=errorbar(common(1).mus(:,1),common(1).mus(:,2),common(1).sigma(:,2)./sqrt(sum(common(1).origbincounts(:,common(1).rep2ind),2)),'.');
e.Color=[0.6 0.6 0.6];
e.MarkerSize=20;
x=linspace(min(common(1).mus(:,1)),max(common(1).mus(:,1)));
mdl = fitlm(common(1).mus(:,1),common(1).mus(:,2));
rsq=mdl.Rsquared.Adjusted;
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
t=sprintf('R^2 = %0.2f\ny = %0.2fx + %0.2f',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
text(3.7,7,t,'FontSize',28)
[t,s]=findSeqsTTest(common(1),1,1e-5);
seqs=common(1).seqs;
% for i=1:length(common)
%     for j=1:length(common(i).seqs)
%         seqs{end+1}=common(i).seqs{j};
%     end
% end
herrorbar(common(1).mus(s,1),common(1).mus(s,2),common(1).sigma(s,1)./sqrt(sum(common(1).origbincounts(s,common(1).rep1ind),2)),common(1).sigma(s,1)./sqrt(sum(common(1).origbincounts(s,common(1).rep1ind),2)),'bo','blue')
errorbar(common(1).mus(s,1),common(1).mus(s,2),common(1).sigma(s,2)./sqrt(sum(common(1).origbincounts(s,common(1).rep2ind),2)),'bo')
% title('All common sequences')
xlabel('mammalian cells')
ylabel('yeast') 
set(gca,'linewidth',2)
set(gca,'fontsize',28)

%%
f=fopen('seqs_yeastmamm.fasta','w');
fprintf(f,'>seq\n%s\n',common(1).seqs{find(s)});
fclose(f);
data=findMotif(common(1),{'(GCTGTC[A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G])([A|C|T|G]*)([A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G]CTGATGAG[A|C|T|G][A|C|T|G][A|C|T|G])([A|C|T|G]*)([A|C|T|G][A|C|T|G][A|C|T|G]CGAAACAGC)'},{'rbz'});
data=findSeqnum(data);
l1l2seqs=[];
l1l2seqsSIG=[];
l1l2seqsNOTSIG=[];
for i=1:length(data.loop2seqnum)
    if length(data.loop2seqnum{i})==5 & s(i)
        l1l2seqsSIG=[l1l2seqsSIG;data.loop2seqnum{(i)}];
    elseif length(data.loop2seqnum{i})==5 & ~s(i)
        l1l2seqsNOTSIG=[l1l2seqsNOTSIG;data.loop2seqnum{(i)}];
    end
end
[entlo,plo]=findEnt(l1l2seqsSIG,{'significantly different'})
[enthi,phi]=findEnt(l1l2seqsNOTSIG,{'not different'})
[dent,dp]=findDeltaEnt(entlo,enthi,plo,phi,{'yeast diff mamm'})








%% END 






