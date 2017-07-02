addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
g=load('~/Documents/YFS/YFSI_gooddataWstruct.mat');
g=g.gooddata;

%%
nohairpin='^\(*\.+\)*$';

loopstruct=g.loop2struct;

hasnohairpin=[find(~cellfun('isempty',regexp(loopstruct,nohairpin)))];
struct2way=g.ensembstruct(find(hasnohairpin));

seqnohairpin={};
loop1nohairpin={};
loop2nohairpin={};
loop1structnohairpin={};
loop2structnohairpin={};
for i=1:length(hasnohairpin)
    seqnohairpin{end+1}=g.seqs{hasnohairpin(i)};
    loop1nohairpin{end+1}=g.loop1{hasnohairpin(i)};
    loop2nohairpin{end+1}=g.loop2{hasnohairpin(i)};
    loop1structnohairpin{end+1}=g.loop1struct{hasnohairpin(i)};
    loop2structnohairpin{end+1}=g.loop2struct{hasnohairpin(i)};
end
musnohairpin=g.mus(find(hasnohairpin));
length(seqnohairpin)

setfig('nohairpin mus');clf
hist(musnohairpin)

%% pull out seqs with specified loop lengths
l1len=7;
l2len=[30 40 50 60];
endslen=[3 3];
load('MyColormaps','mycmap')
DNA={'A','T','C','G'};

%% go through all sequences and convert to matrix of 1 2 3 4
loop1seq=cell(length(g.seqs),1);
loop2seq=loop1seq;
for i=1:length(g.seqs)
    s=g.loop1{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop1seq{i}=s-'0';
    
    s=g.loop2{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop2seq{i}=s-'0';
    
end

g.loop1seqnum=loop1seq;
g.loop2seqnum=loop2seq;

%% Analyses
% instead of having a percentile requirement could we implement 2d struct
% filter?
% muth=g.mus<prctile(g.mus,100);
structfilter=zeros(1,length(g.seqs));
structfilter(hasnohairpin)=1;

[l1l2seqs,l1l2,ttl]=findSeqs(g,l1len,l2len,endslen,structfilter);

% what is the distribution of mus
findParameters(g.mus(l1l2),l1l2seqs);

% find mutual sequence contribution - pairwise mu
[pairwisemu,pairwisesigma]=findPairwiseMu(g.mus(l1l2),l1l2seqs,[],1);
pairwise50=findPairwiseMu(g.mus(l1l2),l1l2seqs,50,[]);
pairwise25=findPairwiseMu(g.mus(l1l2),l1l2seqs,25,[]);
pairwise10=findPairwiseMu(g.mus(l1l2),l1l2seqs,10,[]);

%%
% find entropy
ent=findEnt(l1l2seqs);

% find mutual information
[mI,pJoint,allI]=findMutualInformation(l1l2seqs);

%% display pairwise mu only for positions that matter
% find a cutoff for mI
setfig('mI cutoff');clf
hist(mI(:),length(mI(:))/2)
load('MyColormaps','mycmap')
mIhigh=mI>0.05;
setfig('high mI');clf
imagesc(mIhigh)
mycmap=[[1.00 1.00 1.00];mycmap(end:-1:1,:)];
colormap(mycmap)


%%
% filter pairwise mu for bases with high MI
pairwise_mIhigh25=pairwise25;

for i=1:length(DNA)
    for j=1:length(DNA)
        pairwise_mIhigh25(i:4:end,j:4:end)=mIhigh.*pairwise_mIhigh25(i:4:end,j:4:end);
    end
end

setfig('pairwise 25th percentile mu');clf
imagesc(pairwise_mIhigh25)
for graphing=1
colormap(mycmap)
c=colorbar;
ylabel(c,'\mu')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('pairwise 25th percentile mu')
L = get(gca,'XLim');
xlabelnames=DNA;
a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
NumTicks = length(a)+1;
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{a{:},''})
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'Yticklabel',{a{:},''})
mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);

    for i=1:length(mymajorgrids)
        line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
        line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
    end
end

%% Let's find the consensus and see if there is any shift in library wide mu
[rs cs] = ind2sub(size(mIhigh),find(mIhigh));
pairwise_consensus=zeros(length(rs),2);
for i=1:length(rs)
    if rs(i)>cs(i)
        x=4*(rs(i)-1)+1;
        y=4*(cs(i)-1)+1;
        [r,c]=ind2sub([4 4],find(pairwise_mIhigh25(x:x+3,y:y+3)==min(min(pairwise_mIhigh25(x:x+3,y:y+3)))));
        pairwise_consensus(i,:)=[r c];
    end
end

consensus=[rs cs pairwise_consensus];
consensus=consensus(find(pairwise_consensus(:,1)>0),:)
cmI=[];
for i=1:length(consensus(:,1))
    cmI(end+1)=mI(consensus(i,1),consensus(i,2));
    
end
consensus=[consensus cmI'];
sum(cmI)

%% 
scores=zeros(length(l1l2),1);
for k=1:length(l1l2)
    for i=1:length(consensus(:,1))
        seq=zeros(1,totallooplength);
        seq(consensus(i,1))=consensus(i,3);
        seq(consensus(i,2))=consensus(i,4);
        if sum(l1l2seqs(k,:)==seq)>1
            scores(k)=scores(k)+consensus(i,5);
        end
    end
end

setfig('high scores');clf
hist(scores(scores>0.4),12)
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)
ylabel('frequency')
xlabel('mutual information')

setfig('high scores mus');clf
hist(g.mus(l1l2(scores>0.1)))
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)

g.loop1(l1l2(scores>0.1))
g.loop2(l1l2(scores>0.1))
%% prepare csv output
totallooplength=ttl;
% Incorporate 2-mer based upon loop-interactions
% run analyzePSC first to get consensus (see bottom of that script)

twomer=zeros(length(l1l2),length(consensus(:,1)));
feat2mer={};
for m=1:length(consensus(:,1))
for md=1:4
for nd=1:4
feat2mer{end+1}=sprintf('%0.0f-%0.0f:%s-%s',consensus(m,1),consensus(m,2),DNA{md},DNA{nd});
end
end
end

for i=1:length(l1l2)
    for m=1:length(consensus(:,1))
        for md=1:4
            for nd=1:4
                if l1l2seqs(i,consensus(m,1))==md && l1l2seqs(i,consensus(m,2))==nd
                    twomer(i,16*(m-1)+4*(md-1)+nd)=1;
                end
            end
        end
    end
end

setfig('two mer');clf
bar(sum(twomer))

bitmat=zeros(length(l1l2),length(l1l2seqs(1,:))*4);

for i=1:4
    bitmat(:,i:4:end)=l1l2seqs==i;
end
bitmat=[bitmat twomer];
feats={};
for i=1:l1len
    for j=1:4
    feats{end+1}=sprintf('l1+%0.0f%s',i,DNA{j});
    end
end


for i=1:endslen(1)
    for j=1:4
    feats{end+1}=sprintf('l2+%0.0f%s',i,DNA{j});
    end
end

for i=(endslen(2)):-1:1
    for j=1:4
    feats{end+1}=sprintf('l2-%0.0f%s',i,DNA{j});
    end
end

feats={feats{:} feat2mer{:}};
feats{end+1}='mu';
fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/N7_end3-3_nohairpin_bitmat.csv','w');
csvtitle=sprintf('%s,',feats{:});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',bitmat(i,:));
    ss=sprintf('%s%0.4f',s,g.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end
fclose(fid);

%% END



