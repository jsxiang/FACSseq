addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq/
clear

d=load('FA.mat');
allrbz=d.FA;
x=linspace(min(min(allrbz.VYBmus)),max(max(allrbz.VYBmus)));
labelnames={'no ligand','with ligand'};
numseqs=10;
allrbz(1).motifname='all FA';
allrbzdata=findSEMbyBootstrp(allrbz,100,1:12);
%%
allrbzdata.ttestH=ttest2(allrbzdata.bootcountsR1',allrbzdata.bootcountsR2',1e-2,'both','equal');

%% troubleshoot
q=44;
setfig('troubleshoot');clf
subplot(2,2,1)
hist(allrbzdata.bootcountsR1(q,:))

subplot(2,2,2)
hist(allrbzdata.bootcountsR2(q,:))

subplot(2,2,3)
bar(allrbzdata.scaledbincounts(q,allrbzdata.rep1ind))
subplot(2,2,4)
bar(allrbzdata.scaledbincounts(q,allrbzdata.rep2ind))


%%
setfig('FA switches');clf
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
%%

%% pull out sequences that contain t3hRECK backbone
motifname={'sTRSV','loopI_19nt','loopI_21nt'};
motif={'(GCTGTCACCGG[A|C|T|G])([A|C|T|G]*)([A|C|T|G]CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
       '(GCTGTCACCGG[A|C|T|G])(GCTTGGTACGTTATATTCA)([A|C|T|G]CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
       '(GCTGTCACCGG[A|C|T|G])(TGCTTGGTACGTTATATTCAG)([A|C|T|G]CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
       };
data=findMotif(allrbzdata(1),motif,motifname);
data=findSeqnum(data);
% data=findFoldchange(data,250);


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
findLengthdependence(data(1),mean(data(1).VYBmus,2))
%% compare constrained loops, also do statistical testing
setfig('compare rbz');clf
hold on
legendnames={};
for i=1:length(data)
    subplot(length(data),2,2*i-1)
    histogram(data(i).VYBmus(:,1),linspace(0,5,100),'normalization','probability','edgecolor','none')
    legendnames{end+1}=data(i).motifname;
    title(data(i).motifname,'interpreter','none');
    xlabel('log10(GFP/mCherry)')
    
    
    subplot(length(data),2,2*i)
    prctile(data(i).foldchange,50)
    histogram(data(i).foldchange,linspace(0,400,100),'edgecolor','none')
    xlabel('fold change')
    title(data(i).motifname,'interpreter','none');
    
end
hold off
%%
setfig('FA again');clf
hold on
title('Folinic acid')

plot(data(1).VYBmus(:,1),data(1).VYBmus(:,2),'.','MarkerSize',15)
% plot(data(3).VYBmus(:,1),data(3).VYBmus(:,2),'.','MarkerSize',15)
dscatter(data(1).VYBmus(:,1),data(1).VYBmus(:,2),'BINS',[50 50])
plot(data(1).VYBmus(data(1).switchhits,1),data(1).VYBmus(data(1).switchhits,2),'ro','MarkerSize',8)
x=linspace(min(data(1).VYBmus(:,1)),max(data(1).VYBmus(:,2)));
plot(x,x,'k:')
f=fopen('FAhits.fasta','w');
fprintf(f,'>seq\n%s\n',data(2).seqs{(data(2).switchhits)},data(3).seqs{data(3).switchhits});
fclose(f);
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('log10(GFP/mCherry) no ligand')
ylabel('log10(GFP/mCherry) with ligand')

%% Analyses
i=3;
l2len=[5];
l1len=[21]
endslen=[11 10];

muthlo=data(i).VYBmus(:,1)<prctile(data(i).VYBmus(:,1),50); % filter
muthhi=data(i).VYBmus(:,1)>prctile(data(i).VYBmus(:,1),60); % filter
% muth=~cellfun('isempty',regexp(gooddata.loop1,'^T[A|T|C|G]+'));
[l1l2seqslo,l1l2lo]=findSeqs(data(i),l1len,l2len,endslen,muthlo');
[l1l2seqshi,l1l2hi]=findSeqs(data(i),l1len,l2len,endslen,muthhi');
% seqlogo(data(i).seqs(l1l2lo))
% % what is the distribution of mus
% findParameters(gooddata.VYBmus(l1l2),l1l2seqs);
% find mutual sequence contribution - pairwise mu
[pairwisemulo,pairwisesigmalo]=findPairwiseMu(data(i).VYBmus(l1l2lo,1),l1l2seqslo,10,1,data(i).motifname);
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
[mIhi,pJointhi,allIhi]=findMutualInformation(l1l2seqshi,data(i).motifname);
[deltaMI,deltapJoint]=findDeltaMI(mIlo,mIhi,pJointlo,pJointhi,data(i).motifname);















