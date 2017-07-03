addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
g=load('YFSI_gooddata_noBin1.mat');
gooddata=g.gooddata;
load('MyColormaps','mycmap')
DNA={'A','T','C','G'};

%% go through all sequences and convert to matrix of 1 2 3 4
loop1seq=cell(length(gooddata.seqs),1);
loop2seq=loop1seq;
for i=1:length(gooddata.seqs)
    s=gooddata.loop1{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop1seq{i}=s-'0';
    
    s=gooddata.loop2{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop2seq{i}=s-'0';
    
end

gooddata.loop1seqnum=loop1seq;
gooddata.loop2seqnum=loop2seq;

% %% pull out seqs with specified loop lengths
% l1len=[15:25 30 40 50 60];
% l2len=8;
% 
% endslen=[2 3];
% totallooplength=sum(endslen)+l2len;
% l1ind=zeros(1,length(gooddata.seqs));   
% for i=1:length(l1len)
%     l1ind=l1ind| gooddata.loop1len==l1len(i);
% end
% 
% l2ind=zeros(1,length(gooddata.seqs));
% for i=1:length(l2len)
%     l2ind=l2ind| gooddata.loop2len==l2len(i);
% end
%     
% % let's apply an enrichment to only look at mu's less than the 25th
% % percentile
% % mu25th=gooddata.mus<prctile(gooddata.mus,25);
% % l1l2=find(l1ind.*l2ind.*mu25th);
% 
% l1l2=find(l1ind.*l2ind);
% l1l2seqs=zeros(length(l1l2),totallooplength);
% for i=1:length(l1l2)
%     l1=gooddata.loop1seqnum{l1l2(i)};
%     l1start=l1(1:endslen(1));
%     l1end=l1((end-endslen(2)+1):end);
%     
%     l1l2seqs(i,1:endslen(1))=l1start;
%     l1l2seqs(i,(endslen(1)+1):(endslen(1)+endslen(2)))=l1end;
%    
%     l1l2seqs(i,(sum(endslen)+1):(sum(endslen)+l2len))=gooddata.loop2seqnum{l1l2(i)};
% end
% 
% 
% 
% %% Analyses
% % what is the distribution of mus
% findParameters(gooddata.mus(l1l2),l1l2seqs);
% 
% % find mutual sequence contribution - pairwise mu
% pairwise25=findPairwiseMu(gooddata.mus(l1l2),l1l2seqs,25,[]);
% pairwise10=findPairwiseMu(gooddata.mus(l1l2),l1l2seqs,10,[]);
% [pairwisemu,pairwisesigma]=findPairwiseMu(gooddata.mus(l1l2),l1l2seqs,[],1);
% 
% setfig('pairwise difference');clf
% imagesc(pairwise25-pairwise10)
% 
% for graphing=1
% colormap(mycmap)
% c=colorbar;
% ylabel(c,'\mu')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('pairwise difference')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% NumTicks = length(a)+1;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
% for i=1:length(mymajorgrids)
%     line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%     line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
% end
% end
% 
% % find entropy
% ent=findEnt(l1l2seqs);
% 
% % find mutual information
% [mI,pJoint,allI]=findMutualInformation(l1l2seqs);
% 
% 
% %% display pairwise mu only for positions that matter
% % find a cutoff for mI
% setfig('mI cutoff');clf
% hist(mI(:),length(mI(:))/2)
% set(gca,'fontsize',13)
% set(gca,'linewidth',1.5)
% ylabel('frequency')
% xlabel('mutual information')
% hold on
% plot(0.05, linspace(0,40),'k:','linewidth',2)
% hold off
% 
% load('MyColormaps','mycmap')
% mIhigh=mI>0.05;
% setfig('high mI');clf
% imagesc(mIhigh)
% mycmap=[[1.00 1.00 1.00];mycmap(end:-1:1,:)];
% colormap(mycmap)
% 
% % filter pairwise mu for bases with high MI
% pairwise_mIhigh25=pairwise25;
% 
% for i=1:length(DNA)
%     for j=1:length(DNA)
%         pairwise_mIhigh25(i:4:end,j:4:end)=mIhigh.*pairwise_mIhigh25(i:4:end,j:4:end);
%     end
% end
% 
% setfig('pairwise 25th percentile mu');clf
% imagesc(pairwise_mIhigh25)
% for graphing=1
% colormap(mycmap)
% c=colorbar;
% ylabel(c,'\mu')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('pairwise 25th percentile mu')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% NumTicks = length(a)+1;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
%     for i=1:length(mymajorgrids)
%         line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%         line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
%     end
% end
% 
% % filter pairwise mu for bases with high MI
% pairwise_mIhigh10=pairwise10;
% 
% for i=1:length(DNA)
%     for j=1:length(DNA)
%         pairwise_mIhigh10(i:4:end,j:4:end)=mIhigh.*pairwise_mIhigh10(i:4:end,j:4:end);
%     end
% end
% 
% setfig('pairwise 10th percentile mu');clf
% imagesc(pairwise_mIhigh10)
% for graphing=1
% ncmap=[mycmap(1,:); repmat(mycmap(2,:),round(length(mycmap(:,1))/1.5),1); mycmap(2:end,:)];
% colormap(ncmap)
% c=colorbar;
% ylabel(c,'\mu')
% set(c,'fontsize',13)
% set(c,'linewidth',1.5)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('pairwise 10th percentile mu')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% NumTicks = length(a)+1;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
% for i=1:length(mymajorgrids)
%     line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%     line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
% end
% end
% 
% setfig('pairwise 10th percentile mu adjusted');clf
% p=log(pairwise_mIhigh10.*(1./(pairwise_mIhigh25-pairwise_mIhigh10).^2));
% imagesc(p)
% for graphing=1
%     load('MyColormaps','mycmap')
%     mycmap=[[1.0 1.0 1.0];mycmap];
% colormap(mycmap)
% c=colorbar;
% ylabel(c,'\mu adjusted')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('pairwise 10th percentile mu, adjusted')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% NumTicks = length(a)+1;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
% for i=1:length(mymajorgrids)
%     line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%     line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
% end
% end
% 
% 
% % filter pairwise mu for bases with high MI
% pairwise_mIhigh_mu=pairwisemu;
% 
% for i=1:length(DNA)
%     for j=1:length(DNA)
%         pairwise_mIhigh_mu(i:4:end,j:4:end)=mIhigh.*pairwise_mIhigh_mu(i:4:end,j:4:end);
%     end
% end
% 
% setfig('pairwise mu gaussian');clf
% imagesc(pairwise_mIhigh_mu)
% for graphing=1
% colormap(mycmap)
% c=colorbar;
% ylabel(c,'\mu')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('pairwise mu gaussian')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% NumTicks = length(a)+1;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
% for i=1:length(mymajorgrids)
%     line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%     line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
% end
% end
% 
% 
% 
% %% filter joint probability 
% 
% % filter pairwise mu for bases with high MI
% allI_mIhigh=allI;
% 
% for i=1:length(DNA)
%     for j=1:length(DNA)
%         allI_mIhigh(i:4:end,j:4:end)=mIhigh.*allI_mIhigh(i:4:end,j:4:end);
%     end
% end
% 
% setfig('Joint probability MI High');clf
% allI_mIhigh(allI_mIhigh == 0) = NaN;
% imagesc(allI_mIhigh)
% for graphing=1
% colormap(mycmap)
% c=colorbar;
% ylabel(c,'information')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('information content')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% NumTicks = length(a)+1;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
%     for i=1:length(mymajorgrids)
%         line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%         line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
%     end
% end
% 
% setfig('2 power Joint probability MI High');clf
% imagesc(2.^(allI_mIhigh))
% for graphing=1
% colormap(mycmap)
% c=colorbar;
% ylabel(c,'power 2 information')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('power 2 information content')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% NumTicks = length(a)+1;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
%     for i=1:length(mymajorgrids)
%         line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%         line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
%     end
% end

% 
% %% when loop 2 is AAAGA 
% setfig('loop 2 == AAAGA');clf
% l2AAAGA=find(l1l2seqs(:,end)==1 &l1l2seqs(:,end-1)==4 &l1l2seqs(:,end-2)==1 &l1l2seqs(:,end-3)==1 &l1l2seqs(:,end-4)==1);
% 
% for i=1:length(l2AAAGA)
% seqs=gooddata.seqs{l1l2(l2AAAGA(i))}
% mus=gooddata.mus(l1l2(l2AAAGA(i)))
% gooddata.loop1{l1l2(l2AAAGA(i))}
% gooddata.loop2{l1l2(l2AAAGA(i))}
% subplot(length(l2AAAGA),2,2*i-1)
% bar(gooddata.bincounts(l1l2(l2AAAGA(i)),1:12))
% subplot(length(l2AAAGA),2,2*i)
% bar(gooddata.bincounts(l1l2(l2AAAGA(i)),13:24))
% 
% end
% 
% %% when the sequence contains the theoaptamer 
% setfig('theo');clf
% 
% theo='ATACCAGCATCGTCTTGATGCCCTTGGAAG';
% theoind=zeros(1,length(gooddata.seqs));
% for i=1:length(gooddata.seqs)
%     if regexp(gooddata.seqs{i},'ATACCAGCATCGTCTTGATGCCCTTGGAAG')
%         theoind(i)=1;
%     end
% end
% theoi=find(theoind>0)
% 
% for i=1:length(theoi)
% seqs=gooddata.seqs{theoi(i)}
% mus=gooddata.mus(theoi(i))
% gooddata.loop1{theoi(i)}
% gooddata.loop2{theoi(i)}
% subplot(length(theoi),2,2*i-1)
% bar(gooddata.bincounts(theoi(i),1:12))
% subplot(length(theoi),2,2*i)
% bar(gooddata.bincounts(theoi(i),13:24))
% 
% end
% 
% %% when the sequence contains the beginning and end of theoaptamer 
% setfig('theo');clf
% 
% theo='ATACCAGCATCGTCTTGATGCCCTTGGAAG';
% theoind=zeros(1,length(gooddata.seqs));
% for i=1:length(gooddata.seqs)
%     if gooddata.loop1len(i)>14 & regexp(gooddata.loop1{i},'^AT.*AAG$')
%         theoind(i)=1;
%     elseif gooddata.loop2len(i)>14 & regexp(gooddata.loop2{i},'^AT.*AAG$')
%         theoind(i)=1;
%     end
% end
% theoi=find(theoind>0)
% theostem.mus=[];
% theostem.seqs={};
% theostem.l1={};
% theostem.l2={};
% for i=1:length(theoi)
% theostem.seqs{end+1}=gooddata.seqs{theoi(i)};
% theostem.mus(end+1)=gooddata.mus(theoi(i));
% theostem.l1{end+1}=gooddata.loop1{theoi(i)};
% theostem.l2{end+1}=gooddata.loop2{theoi(i)};
% subplot(length(theoi),2,2*i-1)
% bar(gooddata.bincounts(theoi(i),1:12))
% subplot(length(theoi),2,2*i)
% bar(gooddata.bincounts(theoi(i),13:24))
% 
% end


%% pull out seqs with specified loop lengths
l1len=8;
l2len=[15:25 30 40 50 60];

endslen=[7 8];
totallooplength=l1len+sum(endslen);
l1ind=zeros(1,length(gooddata.seqs));
for i=1:length(l1len)
    l1ind=l1ind| gooddata.loop1len==l1len(i);
end

l2ind=zeros(1,length(gooddata.seqs));
for i=1:length(l2len)
    l2ind=l2ind| gooddata.loop2len==l2len(i);
end


% let's apply an enrichment to only look at mu's less than the 25th
% percentile
mu25th=gooddata.mus<prctile(gooddata.mus,25);
l1l2=find(l1ind.*l2ind.*mu25th);
l1l2all=find(l1ind.*l2ind);
l1l2=find(l1ind.*l2ind);
l1l2seqs=zeros(length(l1l2),totallooplength);
for i=1:length(l1l2)
    l2=gooddata.loop2seqnum{l1l2(i)};
    l2start=l2(1:endslen(1));
    l2end=l2((end-endslen(2)+1):end);
    
    l1l2seqs(i,(l1len+1):(l1len+endslen(1)))=l2start;
    l1l2seqs(i,(l1len+endslen(1)+1):(l1len+endslen(1)+endslen(2)))=l2end;
   
    l1l2seqs(i,1:l1len)=gooddata.loop1seqnum{l1l2(i)};
end
l1l2seqsall=zeros(length(l1l2all),totallooplength);
for i=1:length(l1l2all)
    l2=gooddata.loop2seqnum{l1l2all(i)};
    l2start=l2(1:endslen(1));
    l2end=l2((end-endslen(2)+1):end);
    
    l1l2seqsall(i,(l1len+1):(l1len+endslen(1)))=l2start;
    l1l2seqsall(i,(l1len+endslen(1)+1):(l1len+endslen(1)+endslen(2)))=l2end;
   
    l1l2seqsall(i,1:l1len)=gooddata.loop1seqnum{l1l2all(i)};
end


%% Analyses
% what is the distribution of mus
findParameters(gooddata.mus(l1l2),l1l2seqs);


gooddata.normmus=(gooddata.VYBmus-min(gooddata.VYBmus))./(max(gooddata.VYBmus)-min(gooddata.VYBmus));
% find mutual sequence contribution - pairwise mu
[pairwisemu,pairwisesigma]=findPairwiseMu(gooddata.normmus(l1l2),l1l2seqs,[],1);
pairwise50=findPairwiseMu(gooddata.normmus(l1l2),l1l2seqs,50,[]);
pairwise25=findPairwiseMu(gooddata.normmus(l1l2),l1l2seqs,25,[]);
pairwise10=findPairwiseMu(gooddata.normmus(l1l2),l1l2seqs,10,[]);

setfig('pairwise difference');clf
imagesc(pairwise25-pairwise10)

for graphing=1
colormap(mycmap)
c=colorbar;
ylabel(c,'\mu')
set(c,'fontsize',12)
set(c,'linewidth',1.5)
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
title('pairwise difference')
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

% find entropy
ent=findEnt(l1l2seqs);

% find mutual information
[mI,pJoint,allI]=findMutualInformation(l1l2seqs);


%% display pairwise mu only for positions that matter
% find a cutoff for mI
setfig('mI cutoff');clf
hist(mI(:),length(mI(:))/2)
load('MyColormaps','mycmap')
mIhigh=mI>0.045;
setfig('high mI');clf
imagesc(mIhigh)
mycmap=[[1.00 1.00 1.00];mycmap(end:-1:1,:)];
colormap(mycmap)

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
colormap(mycmap)
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
colormap(mycmap)
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
colormap(mycmap)
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
hist(scores(scores>0.01),12)
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)
ylabel('frequency')
xlabel('mutual information')

setfig('high scores mus');clf
hist(gooddata.mus(l1l2(scores>0.1)))
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)

gooddata.loop1(l1l2(scores>0.1))
gooddata.loop2(l1l2(scores>0.1))


%% convert to bit 
bitmat=zeros(length(l1l2),length(l1l2seqs(1,:))*4);

for i=1:4
    bitmat(:,i:4:end)=l1l2seqs==i;
end

%% Incorporate 2-mer based upon loop-interactions
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
%%
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

bitmat=[bitmat twomer];

fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat_largeII.csv','w');
csvtitle=sprintf('%s,',feats{:});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',bitmat(i,1:end));
    ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end
fclose(fid);

fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat_largeII_2meronly.csv','w');
csvtitle=sprintf('%s,',feats{end-length(feat2mer):end});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',twomer(i,:));
    ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end
fclose(fid);



%% END








