addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq/
clear

d=load('theo.mat');
allrbz=d.theo;
x=linspace(min(min(allrbz.VYBmus)),max(max(allrbz.VYBmus)));
allrbz(1).motifname='all theo';
allrbzdata=findSEMbyBootstrp(allrbz,100,1:12);

allrbzdata.ttestH=ttest2(allrbzdata.bootcountsR1',allrbzdata.bootcountsR2',1e-2,'both','equal');

%%
setfig('Theo switches');clf
% dscatter(allrbzdata.VYBmus(:,1),allrbzdata.VYBmus(:,2),'bins',[20 20])
plot(allrbzdata.VYBmus(:,1),allrbzdata.VYBmus(:,2),'.','Markersize',12)
hold on
h=herrorbar(allrbzdata.VYBmus(:,1),allrbzdata.VYBmus(:,2),allrbzdata.semR1,allrbzdata.semR1,'.',[0.6 0.6 0.6]);
e=errorbar(allrbzdata.VYBmus(:,1),allrbzdata.VYBmus(:,2),allrbzdata.semR2,allrbzdata.semR2,'.','color',[0.6 0.6 0.6]);
e.Color=[0.6 0.6 0.6];
e.MarkerSize=15;
h=herrorbar(allrbzdata.VYBmus(allrbzdata.ttestH==1,1),allrbzdata.VYBmus(allrbzdata.ttestH==1,2),allrbzdata.semR1(allrbzdata.ttestH==1),allrbzdata.semR1(allrbzdata.ttestH==1),'.','blue');
e=errorbar(allrbzdata.VYBmus(allrbzdata.ttestH==1,1),allrbzdata.VYBmus(allrbzdata.ttestH==1,2),allrbzdata.semR2(allrbzdata.ttestH==1),allrbzdata.semR2(allrbzdata.ttestH==1),'b.');

% plot(allrbzdata.VYBmus(allrbzdata.ttestH==1,1),allrbzdata.VYBmus(allrbzdata.ttestH==1,2),'ro','MarkerSize',15)
plot(x,x,'k:')
hold off

%% troubleshoot
q=1;
setfig('troubleshoot');clf
subplot(2,2,1)
hist(allrbzdata.bootcountsR1(q,:))
% xlim([8 11])

subplot(2,2,2)
hist(allrbzdata.bootcountsR2(q,:))
% xlim([8 11])

subplot(2,2,3)
bar(allrbzdata.scaledbincounts(q,allrbzdata.rep1ind))
xlim([4 12.5])

subplot(2,2,4)
bar(allrbzdata.scaledbincounts(q,allrbzdata.rep2ind))
xlim([4 12.5])



%% pull out sequences that contain t3hRECK backbone
motifname={'all ribozymes','sTRSV','t3hRECK'};
motif={'(GCTGTC)([A|C|T|G]*)(CTGATGA)([A|C|T|G]*)(GAAACAGC)',...
       '(GCTGTCACCGG[A|C|T|G])([A|C|T|G]*)([A|C|T|G]CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
       '(GCTGTCCTGCAG)([A|C|T|G]*)(CTGCAGCTGATGAGCTC)([A|C|T|G]*)(GAGCGAAACAGC)',...
       };
data=findMotif(allrbz(1),motif,motifname);
data=findSeqnum(data);
% data=findFoldchange(data,600);


for i=1:length(data)
    switchdata=data(i);
    switchdata.mus=[data(i).mus(:,1) data(i).mus(:,2)];
%     switchdata.sigma=[sqrt(data(i).sigma(:,1).^2+data(i).sigma(:,2).^2) data(i).sigma(:,3)];
    switchdata.sigma=[data(i).sigma(:,1) data(i).sigma(:,2)];
    [t,s]=findSeqsTTest(switchdata,0,1-(1-0.05).^(1/length(switchdata.mus(:,1))),1);
    data(i).ttest=s;
    [f,sh]=findFoldchange(data(i).VYBmus(:,2),data(i).VYBmus(:,1),1);
%     [f,sh]=findFoldchange(data(i).VYBmus(:,2),data(i).combinedVYBmus,2);
    data(i).foldchange=f;
    data(i).switchhits=sh & data(i).ttest';%     data.switchhits=switchhits;

end


%%
for i=1:length(data)
    data(i).combinedVYBmus=mean(data(i).VYBmus,2);
end
%% compare constrained loops, also do statistical testing
setfig('compare rbz');clf
hold on
legendnames={};
for i=1:length(data)
    subplot(length(data),1,i)
    histogram(data(i).mus(:,1),linspace(4,12,500),'normalization','probability','edgecolor','none')
    legendnames{end+1}=data(i).motifname;
    title(data(i).motifname);
    ylim([0,0.1])
end
hold off

% legend(legendnames,'location','best')
%%
setfig('theo again');clf
hold on
title('Theophylline switches')
% plot(data(2).VYBmus(:,1),data(2).VYBmus(:,2),'.','MarkerSize',15)
% plot(data(3).VYBmus(:,1),data(3).VYBmus(:,2),'.','MarkerSize',15)
dscatter(data(1).VYBmus(:,1),data(1).VYBmus(:,2),'BINS',[50 50])
plot(data(1).VYBmus(data(1).switchhits,1),data(1).VYBmus(data(1).switchhits,2),'ro','MarkerSize',8)
x=linspace(min(data(1).VYBmus(:,1)),max(data(1).VYBmus(:,2)));
plot(x,x,'k:')
f=fopen('theohits.fasta','w');
fprintf(f,'>seq\n%s\n',data(2).seqs{(data(2).switchhits)},data(3).seqs{data(3).switchhits});
fclose(f);
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('log10(GFP/mCherry) no ligand')
ylabel('log10(GFP/mCherry) with ligand')
%% 
r=load('rbz.mat');
rbz=r.rbz;
rbzdata=findMotif(rbz(1),motif,motifname);
rbzdata=findSeqnum(rbzdata);

for i=1:length(rbzdata)
    switchdata=rbzdata(i);
    switchdata.mus=[rbzdata(i).mus(:,1) rbzdata(i).mus(:,2)];
%     switchdata.sigma=[sqrt(data(i).sigma(:,1).^2+data(i).sigma(:,2).^2) data(i).sigma(:,3)];
    switchdata.sigma=[rbzdata(i).sigma(:,1) rbzdata(i).sigma(:,2)];
    [t,s]=findSeqsTTest(switchdata,0,1-(1-0.05).^(1/length(switchdata.mus(:,1))),1);
    rbzdata(i).ttest=s;
    [f,sh]=findFoldchange(rbzdata(i).VYBmus(:,2),rbzdata(i).VYBmus(:,1),1);
%     [f,sh]=findFoldchange(rbzdata(i).VYBmus(:,3),rbzdata(i).combinedVYBmus,2);
    rbzdata(i).foldchange=f;
    rbzdata(i).switchhits=sh & rbzdata(i).ttest';%     data.switchhits=switchhits;
end

%%


sTRSV='GCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGC';
g814='GCTGTCACCGGATGTTTCCGGTCTGATGAGTCCATAAGGACGAAACAGC';
g833='GCTGTCACCGGATACTTCCGGTCTGATGAGTCCCAGAGGACGAAACAGC';
g862='GCTGTCACCGGATGCATCCGGTCTGATGAGTCCCGCGGGACGAAACAGC';
L2b9a1='GCTGTCACCGGAATCAAGGTCCGGTCTGATGAGTCCGTTGTCCAATACCAGCATCGTCTTGATGCCCTTGGCAGTGGATGGGGACGGAGGACGAAACAGC';
Rs2='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAAAGGACGAAACAGC';
Theo1041='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAGAGGACGAAACAGC';
sTRSVctl='GCTGTCACCGGATGTGCTTTCCGGTACGTGAGGTCCGTGAGGACGAAACAGC';

theolibcontrols={sTRSV,g814,g833,L2b9a1,Theo1041,Rs2};
ctrlnames={'sTRSV','g814','g833','L2b9a1','Theo1041','Rs2'};
for k=1:length(theolibcontrols)
    for j=1
    s=regexp(data(1).seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    ctrlmus=data(1).VYBmus(ctrl,:);
    plot(ctrlmus(:,1),ctrlmus(:,2),'+','linewidth',2,'MarkerSize',12)
    s=regexp(rbzdata(1).seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    ctrlmus=rbzdata(1).VYBmus(ctrl,:);
    plot(ctrlmus(:,1),ctrlmus(:,2),'+','linewidth',2,'MarkerSize',12)
    end
end
legend({'library','hits','1:1',ctrlnames{:}},'location','best')
%% ttest
setfig('ttest');clf
ttestmat=[];
for i=1:length(data)
    for j=1:length(data)
        subplot(length(data),length(data),length(data)*(i-1)+j)
        ml=min(length(data(j).mus),length(data(i).mus));
        inds1=randperm(ml);
        inds2=randperm(ml);
        if i==j
            histogram(mean(data(i).mus(inds1,:),2)-mean(data(j).mus(inds2,:),2),'normalization','probability','edgecolor','none','facecolor',[0.6 0.6 0.6]);
        else
            histogram(mean(data(i).mus(inds1,:),2)-mean(data(j).mus(inds2,:),2),'normalization','probability','edgecolor','none');
        end
        xlim([-8 8])
        [h,p]=ttest(data(i).mus(1:ml,1),data(j).mus(1:ml,1));
        ttestmat(i,j)=p;
        if i==length(data)
            xlabel(data(j).motifname)
            if j==1
               ylabel(data(i).motifname)
            end
        elseif j==1
            ylabel(data(i).motifname)
        end
    end
end

setfig('heat');clf
imagesc(log10(ttestmat))
colorbar

%%
findLengthdependence(data(1),mean(data(1).mus,2))
%%
%% Analyses
%% Analyses
i=2;
l2len=[30];
l1len=[5]
endslen=[15 15];

muthlo=data(i).combinedVYBmus<prctile(data(i).combinedVYBmus,25); % filter
muthhi=data(i).combinedVYBmus>prctile(data(i).combinedVYBmus,60); % filter
% muth=~cellfun('isempty',regexp(gooddata.loop1,'^T[A|T|C|G]+'));
[l1l2seqslo,l1l2lo]=findSeqs(data(i),l1len,l2len,endslen,muthlo');
[l1l2seqshi,l1l2hi]=findSeqs(data(i),l1len,l2len,endslen,muthhi');
% seqlogo(data(i).seqs(l1l2lo))
% % what is the distribution of mus
% findParameters(gooddata.VYBmus(l1l2),l1l2seqs);
% find mutual sequence contribution - pairwise mu
[pairwisemulo,pairwisesigmalo]=findPairwiseMu(data(i).combinedVYBmus(l1l2lo),l1l2seqslo,10,1,data(i).motifname);
% pairwise50=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,50,[]);
% pairwise25=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,25,[]);
% pairwise10=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,10,[]);

set(gca,'fontsize',14)

% find entropy
[entlo,plo]=findEnt(l1l2seqslo,data(i).motifname);
[enthi,phi]=findEnt(l1l2seqshi,data(i).motifname);

[dent,dp]=findDeltaEnt(entlo,enthi,plo,phi,data(i).motifname);

% find mutual information
[mIlo,pJointlo,allIlo]=findMutualInformation(l1l2seqslo,data(i).motifname);
% [mIhi,pJointhi,allIhi]=findMutualInformation(l1l2seqshi,data(i).motifname);
[deltaMI,deltapJoint]=findDeltaMI(mIlo,mIhi,pJointlo,pJointhi,data(i).motifname);








