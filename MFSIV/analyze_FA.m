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
%% pull out sequences that contain t3hRECK backbone
motifname={'sTRSV','loopI_19nt_4bp','loopI_19nt_3bp','loopI_21nt'};
motif={'(GCTGTCACCGG[A|C|T|G])([A|C|T|G]*)([A|C|T|G]CCGGTCTGATGAGTC)([A|C|T|G]*)(GACGAAACAGC)',...
       '(GCTGTCACCGG[A|C|T|G])(GCTTGGTACGTTATATTCA)([A|C|T|G]CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
       '(GCTGTCACCGG[A|C|T|G])(GCTTGGTACGTTATATTCA)([A|C|T|G]CCGGTCTGATGAGTC)([A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G])(GACGAAACAGC)',...
       '(GCTGTCACCGG[A|C|T|G])(TGCTTGGTACGTTATATTCAG)([A|C|T|G]CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
       };
data=findMotif(allrbz(1),motif,motifname);
data=findSeqnum(data);

for i=1:length(data)
    data(i).combinedVYBmus=mean(data(i).VYBmus(:,[1 3]),2);
    switchdata=data(i);
    switchdata.rep2ind=data(i).rep4ind;
%     switchdata.mus=[mean(data(i).mus(:,[1 2]),2) data(i).mus(:,4)];
    switchdata.mus=[data(i).mus(:,1) data(i).mus(:,4)];
%     switchdata.sigma=[sqrt(data(i).sigma(:,1).^2+data(i).sigma(:,2).^2) data(i).sigma(:,3)];
    switchdata.sigma=[data(i).sigma(:,1) data(i).sigma(:,4)];
    [t,s]=findSeqsTTest(switchdata,0,1-(1-0.05).^(1/length(switchdata.mus(:,1))),1);
    data(i).ttest=s;
%     data(i).switchhits=findFoldchange(data(i).VYBmus(:,4),data(i).combinedVYBmus,1) & data(i).ttest';
    [f,sh]=findFoldchange(data(i).VYBmus(:,4),data(i).combinedVYBmus,2.5);
    data(i).foldchange=f;
    data(i).switchhits=sh;
%     data.switchhits=switchhits;

end
% data=findFoldchange(data,1);
% for i=1:length(data)
%     data(i).switchhits=data(i).switchhits & data(i).ttest';
% end
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
% plot(data(2).combinedVYBmus(:,1),data(2).VYBmus(:,4),'.','MarkerSize',15)
% plot(data(3).combinedVYBmus(:,1),data(3).VYBmus(:,4),'.','MarkerSize',15)
dscatter(data(1).combinedVYBmus(:,1),data(1).VYBmus(:,4),'BINS',[50 50])
plot(data(1).combinedVYBmus(data(1).switchhits),data(1).VYBmus(data(1).switchhits,4),'bo','MarkerSize',8)
plot(data(1).combinedVYBmus(data(1).ttest),data(1).VYBmus(data(1).ttest,4),'ro','MarkerSize',8)
x=linspace(min(data(1).VYBmus(:,1)),max(data(1).VYBmus(:,2)));
plot(x,x,'k:')
f=fopen('FAhits.fasta','w');
fprintf(f,'>seq\n%s\n',data(2).seqs{(data(2).switchhits)},data(3).seqs{data(3).switchhits});
fclose(f);
set(gca,'linewidth',1.5)
set(gca,'fontsize',16)
title('Folinic acid switches')
xlabel('log10(GFP/mCherry) no ligand')
ylabel('log10(GFP/mCherry) with ligand')

%% Analyses
i=2;
l2len=[5];
l1len=[19];
endslen=[10 9];

muthlo=data(i).VYBmus(:,1)<prctile(data(i).VYBmus(:,1),20); % filter
muthhi=data(i).VYBmus(:,1)>prctile(data(i).VYBmus(:,1),80); % filter
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

[ee,ep]=findEntEnrich(entlo,enthi,plo,phi,data(i).motifname);

% find mutual information
[mIlo,pJointlo,allIlo]=findMutualInformation(l1l2seqslo,data(i).motifname);
[mIhi,pJointhi,allIhi]=findMutualInformation(l1l2seqshi,data(i).motifname);
[deltaMI,deltapJoint]=findDeltaMI(mIlo,mIhi,pJointlo,pJointhi,data(i).motifname);
%%
setfig('loopI_19nt_4bp: \DeltaMI')
axis([21 26 21 26]-1.5)
set(gca,'YTicklabels',20:24)
set(gca,'YTick',20:24)

setfig('loopI_19nt_4bp: Entropy enrichment')
set(gca,'XTick',20:24)
xlim([19.5 24.5])

setfig('loopI_19nt_4bp: pairwise mu')
load('MyColormaps','mycmap')
mycmap=mycmap(end:-1:1,:);

mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
axis([19*4 (19*4)+20 19*4 (19*4)+20]+0.5)
set(gca,'YTick',(19*4):((19*4)+20))
set(gca,'YTickLabels',repmat({'G','A','U','C'},1,5))
set(gca,'XTick',(19*4):((19*4)+20))
set(gca,'XTickLabels',repmat({'G','A','U','C'},1,5))



%%
setfig('loopI_21nt: \DeltaMI')
axis([21 26 21 26]+0.5)
set(gca,'YTicklabels',22:26)
set(gca,'YTick',22:26)

setfig('loopI_21nt: Entropy enrichment')
set(gca,'XTick',22:26)
xlim([21.5 26.5])

setfig('loopI_21nt: pairwise mu')
load('MyColormaps','mycmap')
mycmap=mycmap(end:-1:1,:);

mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
axis([21*4 (21*4)+20 21*4 (21*4)+20]+0.5)
set(gca,'YTick',(21*4):((21*4)+20))
set(gca,'YTickLabels',repmat({'G','A','U','C'},1,5))
set(gca,'XTick',(21*4):((21*4)+20))
set(gca,'XTickLabels',repmat({'G','A','U','C'},1,5))



%% where are the switches

setfig('switches');clf

for i=2:length(data)
    subplot(3,3,3*(i-2)+1)
    title('\mu-')
    hold on
    h=histogram(data(i).VYBmus(:,1));
    h.Normalization='probability';
    xlim([-1.5 1.5])
    
    subplot(3,3,3*(i-2)+2)
    title('\mu+')
    hold on
    h=histogram(data(i).VYBmus(:,4));
    h.Normalization='probability';
    xlim([-1.5 1.5])

    subplot(3,3,3*(i-2)+3)
    title('\Delta\mu')
    hold on
    h=histogram(data(i).VYBmus(:,4)-data(i).VYBmus(:,1));
    h.Normalization='probability';
    xlim([-1.5 1.5])
end

%% 

setfig('fold change');clf
hold on
histogram(data(2).VYBmus(:,4)-data(2).VYBmus(:,1));
histogram(data(4).VYBmus(:,4)-data(4).VYBmus(:,1));











