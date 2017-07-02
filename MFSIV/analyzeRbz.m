addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq/
clear

d=load('rbz.mat');
allrbz=d.rbz;
x=linspace(min(min(allrbz.VYBmus)),max(max(allrbz.VYBmus)));
labelnames={'1','2'};
numseqs=10;
allrbz(1).motifname='all ribozymes';
%% pull out sequences that contain t3hRECK backbone
motifname={'sTRSV','t3hRECK','unna','unna-t3hRECK','ungna','ungna-t3hRECK','ungcna','ungcna-t3hRECK','no tert motif','unna\{ungna}','ungna\{ungcna}'};
motif={'(GCTGTCACCGGA)([A|C|T|G]*)(TCCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)',...
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
    data(i).combinedVYBmus=mean(data(i).VYBmus(:,[1 3]),2);
end
%% compare constrained loops, also do statistical testing
setfig('compare rbz');clf
hold on
legendnames={};
for i=1:length(data)
    subplot(length(data),1,i)
    histogram(mean(data(i).VYBmus,2),linspace(-1,5,500),'normalization','probability','edgecolor','none')
    legendnames{end+1}=data(i).motifname;
    title(data(i).motifname);
%     ylim([0,0.1])
end
hold off

% legend(legendnames,'location','best')

%% ttest
setfig('ttest');clf
ttestmat=[];
for i=1:length(data)
    for j=1:length(data)
        subplot(length(data),length(data),length(data)*(i-1)+j)
        ml=min(length(data(j).VYBmus),length(data(i).VYBmus));
        inds1=randperm(ml);
        inds2=randperm(ml);
        if i==j
            histogram(mean(data(i).VYBmus(inds1,:),2)-mean(data(j).VYBmus(inds2,:),2),'normalization','probability','edgecolor','none','facecolor',[0.6 0.6 0.6]);
        else
            histogram(mean(data(i).VYBmus(inds1,:),2)-mean(data(j).VYBmus(inds2,:),2),'normalization','probability','edgecolor','none');
        end
        xlim([-8 8])
        [h,p]=ttest(data(i).VYBmus(1:ml,1),data(j).VYBmus(1:ml,1));
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
findLengthdependence(data(1),mean(data(1).VYBmus,2))

%% save a rbz file

% ttestsig=find(~data(i).ttestdiff(l1l2lo,2));
seqs2save={};
for k=1:length(l1l2lo)

    seqs2save{end+1}=regexprep(data(i).seqs{l1l2lo(k)},'T','U');

end

f=fopen('allsigrbz.fasta','w');
fprintf(f,'>seq\n%s\n',seqs2save{:});
fclose(f);



%% Analyses
i=1;
l1len=[6 7 8];
l2len=[30 30];
endslen=[3 3];

muthlo=data(i).combinedVYBmus<prctile(data(i).combinedVYBmus,20); % filter
muthhi=data(i).combinedVYBmus>prctile(data(i).combinedVYBmus,80); % filter
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
[eent,ep]=findEntEnrich(entlo,enthi,plo,phi,data(i).motifname);


% find mutual information
[mIlo,pJointlo,allIlo]=findMutualInformation(l1l2seqslo,data(i).motifname);
[mIhi,pJointhi,allIhi]=findMutualInformation(l1l2seqshi,data(i).motifname);
[deltaMI,deltapJoint]=findDeltaMI(mIlo,mIhi,pJointlo,pJointhi,data(i).motifname);

%%
setfig('sTRSV: \DeltaMI')
hold on
plot(ones(100,1)*6.5,linspace(0.5,12.5)',':','color',[0.5 0.5 0.5],'linewidth',2)
plot(linspace(0.5,12.5)',ones(100,1)*6.5,':','color',[0.5 0.5 0.5],'linewidth',2)
set(gca,'YTick',1:12)
set(gca,'YTickLabels',{'I1','I2','I3','I-3','I-2','I-1', 'II1','II2','II3','II-3','II-2','II-1'})

set(gca,'XTick',1:12)
set(gca,'XTickLabels',{'I1','I2','I3','I-3','I-2','I-1', 'II1','II2','II3','II-3','II-2','II-1'})

setfig('sTRSV: Entropy enrichment')
hold on
plot(ones(100,1)*6.5,linspace(-0.03,0.03)',':','color',[0.5 0.5 0.5],'linewidth',2)
set(gca,'XTick',1:12)


setfig('sTRSV: pairwise mu')
load('MyColormaps','mycmap')
mycmap=mycmap(end:-1:1,:);

mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
set(gca,'YTick',1:(12*4))
set(gca,'YTickLabels',repmat({'A','U','C','G'},1,5))
set(gca,'XTick',1:(12*4))
set(gca,'XTickLabels',repmat({'A','U','C','G'},1,5))



