addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq/
clear

d=load('theo.mat');
allrbz=d.theo;
x=linspace(min(min(allrbz.mus)),max(max(allrbz.mus)));
labelnames={'1','2'};
numseqs=10;
allrbz(1).motifname='all theo';
%% pull out sequences that contain t3hRECK backbone
motifname={'sTRSV','t3hRECK','unna','unna-t3hRECK','ungna','ungna-t3hRECK','ungcna','ungcna-t3hRECK','no tert motif','unna\{ungna}','ungna\{ungcna}'};
motif={'(GCTGTCACCGG[A|C|T|G])([A|C|T|G]*)([A|C|T|G]CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
       '(GCTGTCCTGCAG)([A|C|T|G]*)(CTGCAGCTGATGAGCTC)([A|C|T|G]*)(GAGCGAAACAGC)',...
       '(GCTGTCACCGGA)(T[A|C|T|G]*)(TCCGGTCTGATGAGTCC)([A|C|T|G]*A)(GGACGAAACAGC)',...
       '(GCTGTCCTGCAG)(T[A|C|T|G]*)(CTGCAGCTGATGAGCTC)([A|C|T|G]*A)(GAGCGAAACAGC)',...
       '(GCTGTCACCGGA)(T[A|C|T|G]*G)(TCCGGTCTGATGAGTCC)([A|C|T|G]*A)(GGACGAAACAGC)',...
       '(GCTGTCCTGCAG)(T[A|C|T|G]*G)(CTGCAGCTGATGAGCTC)([A|C|T|G]*A)(GAGCGAAACAGC)',...
       '(GCTGTCACCGGA)(T[A|C|T|G]*G)(TCCGGTCTGATGAGTCC)(C[A|C|T|G]*A)(GGACGAAACAGC)',...
       '(GCTGTCCTGCAG)(T[A|C|T|G]*G)(CTGCAGCTGATGAGCTC)(C[A|C|T|G]*A)(GAGCGAAACAGC)',...
       '(GCTGTCACCGGA)([A|C|G][A|C|T|G]*[A|T|C])(TCCGGTCTGATGAGTCC)([A|T|G][A|C|T|G]*[T|C|G])(GGACGAAACAGC)',...
       '(GCTGTCACCGGA)(T[A|C|T|G]*[A|T|C])(TCCGGTCTGATGAGTCC)([A|C|T|G]*A)(GGACGAAACAGC)',...
       '(GCTGTCACCGGA)(T[A|C|T|G]*[A|T|C])(TCCGGTCTGATGAGTCC)([A|T|G][A|C|T|G]*A)(GGACGAAACAGC)'};
data=findMotif(allrbz(1),motif,motifname);
data=findSeqnum(data);
%%
for i=1:length(data)
    data(i).combinedmus=mean(data(i).mus,2);
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
i=1;
l2len=[5];
l1len=[30]
endslen=[15 15];

muthlo=data(i).mus(:,1)<prctile(data(i).mus(:,1),25); % filter
muthhi=data(i).mus(:,1)>prctile(data(i).mus(:,1),75); % filter
% muth=~cellfun('isempty',regexp(gooddata.loop1,'^T[A|T|C|G]+'));
[l1l2seqslo,l1l2lo]=findSeqs(data(i),l1len,l2len,endslen,muthlo');
[l1l2seqshi,l1l2hi]=findSeqs(data(i),l1len,l2len,endslen,muthhi');
seqlogo(data(i).seqs(l1l2lo))
% % what is the distribution of mus
% findParameters(gooddata.VYBmus(l1l2),l1l2seqs);
% find mutual sequence contribution - pairwise mu
[pairwisemulo,pairwisesigmalo]=findPairwiseMu(data(i).combinedmus(l1l2lo),l1l2seqslo,10,1,data(i).motifname);
% pairwise50=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,50,[]);
% pairwise25=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,25,[]);
% pairwise10=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,10,[]);

set(gca,'fontsize',14)

% find entropy
[entlo,plo]=findEnt(l1l2seqslo,data(i).motifname);
[enthi,phi]=findEnt(l1l2seqshi,data(i).motifname);
[dent,dp]=findDeltaEnt(entlo,enthi,plo,phi,data(i).motifname);
% find mutual information
[mI,pJoint,allI]=findMutualInformation(l1l2seqslo,data(i).motifname);













