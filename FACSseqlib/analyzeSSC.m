addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
g=load('YFSI_gooddataWstruct.mat');
g=g.gooddata;
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


%% pull out seqs with specified loop lengths
l1len=7;
% l2len=[15:25 30 40 50 60];
l2len=[15];
endslen=[3 3];

%% Analyses
% instead of having a percentile requirement could we implement 2d struct
% filter?
% muth=g.mus<prctile(g.mus,100);
structfilter=zeros(1,length(g.seqs));
structfilter(has2way(symb(sym2bulge)))=1;

[l1l2seqs,l1l2,ttl]=findSeqs(g,l1len,l2len,endslen,structfilter);

% what is the distribution of mus
findParameters(g.mus(l1l2),l1l2seqs);

% find mutual sequence contribution - pairwise mu
[pairwisemu,pairwisesigma]=findPairwiseMu(g.mus(l1l2),l1l2seqs,[],1);
pairwise50=findPairwiseMu(g.mus(l1l2),l1l2seqs,50,[]);
pairwise25=findPairwiseMu(g.mus(l1l2),l1l2seqs,25,[]);
pairwise10=findPairwiseMu(g.mus(l1l2),l1l2seqs,10,[]);

%% prepare csv output
totallooplength=ttl;
%% convert to bit 
%% Incorporate 2-mer based upon loop-interactions
% run analyzePSC first to get consensus (see bottom of that script)
consensus =[
    12     2     1     1
    12     6     2     3
    11     9     1     2];

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
fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/N7_end3-3_2bulgetwowayjunc_bitmat.csv','w');
csvtitle=sprintf('%s,',feats{:});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',bitmat(i,:));
    ss=sprintf('%s%0.4f',s,g.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end
fclose(fid);

% outcsvcell=cell(length(l1l2),totallooplength+1);
% for i=1:length(l1l2)
%     seqindigits=sprintf('%0.0f',l1l2seqs(i,:));
%     s=regexprep(seqindigits,'1','A');
%     s=regexprep(s,'2','T');
%     s=regexprep(s,'3','C');
%     s=regexprep(s,'4','G');
%     for j=1:totallooplength
%         outcsvcell{i,j}=s(j);
%     end
%     outcsvcell{i,totallooplength+1}=sprintf('%0.4f',g.mus(l1l2(i)));
% end
% 
% 
% csvtitle=strcat(sprintf('l1+%0.0f,',1:l1len),...
%                 sprintf('l2+%0.0f,',1:endslen(1)),...
%                 sprintf('l2-%0.0f,',(endslen(2)):-1:1),...
%                 sprintf('mu\n'));
% 
% fid = fopen('~/Documents/YFS/N7_end3-3_2bulgetwowayjunc.csv','w');
% fprintf(fid,'%s\n',csvtitle);
% for i=1:length(l1l2)
%     s=strcat(sprintf('%s,',outcsvcell{i,1:end-1}),outcsvcell(i,end));
%     fprintf(fid,'%s\n',s{1});
% end
% fclose(fid);


%%
% find entropy
ent=findEnt(l1l2seqs);

% find mutual information
[mI,pJoint,allI]=findMutualInformation(l1l2seqs);

%%
muth=(g.mus<prctile(g.mus,25)).*structfilter;
[l1l2seqs_25,l1l2_25]=findSeqs(g,l1len,l2len,endslen,muth);

muth=(g.mus>prctile(g.mus,75)).*structfilter;
[l1l2seqs_75,l1l2_75]=findSeqs(g,l1len,l2len,endslen,muth);

n=mean([length(l1l2_25) length(l1l2_75)])*4

% analyses for the fastest cleaving 25th percentile
% 
findParameters(g.mus(l1l2_25),l1l2seqs_25);

% find entropy
[ent25,pAll25]=findEnt(l1l2seqs_25);
% find mutual information
[mI_25,pJoint_25,allI_25]=findMutualInformation(l1l2seqs_25);

% analyses for the slowest cleaving 25th percentile
findParameters(g.mus(l1l2_75),l1l2seqs_75);

[ent75,pAll75]=findEnt(l1l2seqs_75);
[mI_75,pJoint_75,allI_75]=findMutualInformation(l1l2seqs_75);

% find difference between entropy
deltaEnt=ent25-ent75;
deltaPAll=pAll25-pAll75;
setfig('\Deltaentropy');clf
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
colormap('jet')
legend('A','T','C','G','Location','Best')
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
ylabel('entropy')
xlabel('position')

setfig('\Deltaprobability');clf
imagesc(deltaPAll')
colormap(mycmap)
c=colorbar;
ylabel(c,'\Deltaprobability')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('\Deltaprobability')

L = get(gca,'XLim');
xlabelnames=DNA;
NumTicks = length(xlabelnames)+5;
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{'','A','','T','','C','','G',''})
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)

% difference between max probability and min probability at each position
diffP=max(deltaPAll)-min(deltaPAll);


%% find difference between the MI, pJoint and All_Information

detlaAllI=allI_25-allI_75;
setfig('delta all information');clf
% mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
imagesc(detlaAllI)
c=colorbar;
ylabel(c,'position specific information')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('position specific information')
L = get(gca,'XLim');
xlabelnames=DNA;
a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs_25(1,:)));

NumTicks = length(a)+1;

set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{a{:},''})
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'Yticklabel',{a{:},''})

mymajorgrids=0.5:4:(length(l1l2seqs_25(1,:))*4+0.5);

for i=1:length(mymajorgrids)
    line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
    line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
end

%%
deltaMI=mI_25-mI_75;
setfig('delta MI');clf
% mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
imagesc(deltaMI)
c=colorbar;
ylabel(c,'\DeltaMI')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('\DeltaMI')

%% pairwise delta mu
[pairwisemu25,pairwisesigma25]=findPairwiseMu(g.mus(l1l2_25),l1l2seqs_25,[],1);
[pairwisemu75,pairwisesigma75]=findPairwiseMu(g.mus(l1l2_75),l1l2seqs_75,[],1);

deltapairwisemu=pairwisemu25-pairwisemu75;

setfig('delta pairwise mu information');clf
% mycmap=[[1.0 1.0 1.0];mycmap];
colormap(mycmap)
imagesc(deltapairwisemu)
c=colorbar;
ylabel(c,'\Delta\mu')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('\Deltapairwise\mu')
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

%%
%% display pairwise mu only for positions that matter

% find a cutoff for mI
setfig('mI cutoff');clf
hist(deltaMI(:),length(deltaMI(:))/2)
xlabel('\DeltaMI')
ylabel('Frequency')
load('MyColormaps','mycmap')
mIhigh=deltaMI>0.05;
setfig('high mI');clf
imagesc(mIhigh)
% mycmap=[[1.00 1.00 1.00];mycmap(end:-1:1,:)];
% colormap(mycmap)

% filter pairwise mu for bases with high MI
pairwise_mIhigh25=deltapairwisemu;

for i=1:length(DNA)
    for j=1:length(DNA)
        pairwise_mIhigh25(i:4:end,j:4:end)=mIhigh.*pairwise_mIhigh25(i:4:end,j:4:end);
    end
end

setfig('pairwise 25th percentile mu');clf
imagesc(pairwise_mIhigh25)
for graphing=1
colormap([[1.00 1.00 1.00];mycmap(end:-1:1,:);[1.00 1.00 1.00]])
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

% filter pairwise mu for bases with high MI
pairwise_mIhigh10=pairwise10;

for i=1:length(DNA)
    for j=1:length(DNA)
        pairwise_mIhigh10(i:4:end,j:4:end)=mIhigh.*pairwise_mIhigh10(i:4:end,j:4:end);
    end
end

setfig('pairwise 10th percentile mu');clf
imagesc(pairwise_mIhigh10)
for graphing=1
% colormap(mycmap)
c=colorbar;
ylabel(c,'\mu')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('pairwise 10th percentile mu')
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

setfig('pairwise 10th percentile mu adjusted');clf
p=log(pairwise_mIhigh10.*(1./(pairwise_mIhigh25-pairwise_mIhigh10).^2));
imagesc(p)
for graphing=1
    load('MyColormaps','mycmap')
    mycmap=[[1.0 1.0 1.0];mycmap];
% colormap(mycmap)
c=colorbar;
ylabel(c,'\mu adjusted')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('pairwise 10th percentile mu, adjusted')
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


% filter pairwise mu for bases with high MI
pairwise_mIhigh_mu=pairwisemu;

for i=1:length(DNA)
    for j=1:length(DNA)
        pairwise_mIhigh_mu(i:4:end,j:4:end)=mIhigh.*pairwise_mIhigh_mu(i:4:end,j:4:end);
    end
end

setfig('pairwise mu gaussian');clf
imagesc(pairwise_mIhigh_mu)
for graphing=1
% colormap(mycmap)
c=colorbar;
ylabel(c,'\mu')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('pairwise mu gaussian')
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
    cmI(end+1)=deltaMI(consensus(i,1),consensus(i,2));
    
end
consensus=[consensus cmI'];
sum(cmI)
cmI(isnan(cmI))=0;
%% 
scores=zeros(length(l1l2),1);
for k=1:length(l1l2)
    for i=1:length(consensus(:,1))
        seq=zeros(1,l1len+l2len);
        seq(consensus(i,1))=consensus(i,3);
        seq(consensus(i,2))=consensus(i,4);
        if sum(l1l2seqs(k,:)==seq)>1
            scores(k)=scores(k)+consensus(i,5);
        end
    end
end

setfig('high scores');clf
hist(scores(scores>0.01),12)
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)
ylabel('frequency')
xlabel('mutual information')

setfig('high scores mus');clf
hist(g.mus(l1l2(scores>0.1)))
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)

setfig('mus');clf
hist(g.mus(l1l2))
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)

g.loop1(l1l2(scores>0.1))
g.loop2(l1l2(scores>0.1))

mean(g.mus(l1l2(scores>0.1)))
mean(g.mus(l1l2))

%% 
endisA=l1l2seqs(:,end)==ones(size(l1l2'));
setfig('end is A');clf
hist(g.mus(l1l2(find(endisA))))
mean(g.mus(l1l2(find(endisA))))
%% END








