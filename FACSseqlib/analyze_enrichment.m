addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
% g=load('YFSI_gooddataWstruct.mat');
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


%% pull out seqs with specified loop lengths
% l1len=[5 6 7 8];
% l2len=[15:25 30 40 50 60];
% endslen=[3 3];

l1len=[5];
l2len=[15:25 30 40 50 60];
endslen=[7 8];


%% Analyses

muth=gooddata.VYBmus<prctile(gooddata.VYBmus,100); % filter
% muth=~cellfun('isempty',regexp(gooddata.loop1,'^T[A|T|C|G]+'));
[l1l2seqs,l1l2]=findSeqs(gooddata,l1len,l2len,endslen,muth);

% what is the distribution of mus
findParameters(gooddata.VYBmus(l1l2),l1l2seqs);

% find mutual sequence contribution - pairwise mu
[pairwisemu,pairwisesigma]=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,10,1);
% pairwise50=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,50,[]);
% pairwise25=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,25,[]);
% pairwise10=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,10,[]);

set(gca,'fontsize',14)
%%
% find entropy
[ent,p]=findEnt(l1l2seqs);

% find mutual information
[mI,pJoint,allI]=findMutualInformation(l1l2seqs);


%%
% filter=~cellfun('isempty',regexp(gooddata.loop1,'^T[A|T|C|G]+G$'));
filter=~cellfun('isempty',regexp(gooddata.loop1,'[A|T|C|G]+'));

muth=(gooddata.VYBmus<prctile(gooddata.VYBmus,50)) & filter;
[l1l2seqs_25,l1l2_25]=findSeqs(gooddata,l1len,l2len,endslen,muth);


muth=(gooddata.VYBmus>prctile(gooddata.VYBmus,50)) & filter;
[l1l2seqs_75,l1l2_75]=findSeqs(gooddata,l1len,l2len,endslen,muth);

n=mean([length(l1l2_25) length(l1l2_75)])*4

% analyses for the fastest cleaving 25th percentile
% 
findParameters(gooddata.VYBmus(l1l2_25),l1l2seqs_25);

% find entropy
[ent25,pAll25]=findEnt(l1l2seqs_25);

%%

% find mutual information
[mI_25,pJoint_25,allI_25]=findMutualInformation(l1l2seqs_25);

% analyses for the slowest cleaving 25th percentile
findParameters(gooddata.VYBmus(l1l2_75),l1l2seqs_75);

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
[pairwisemu25,pairwisesigma25]=findPairwiseMu(gooddata.VYBmus(l1l2_25),l1l2seqs_25,10,1);

%%
[pairwisemu50,pairwisesigma50]=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,10,1);


%%
[pairwisemu75,pairwisesigma75]=findPairwiseMu(gooddata.VYBmus(l1l2_75),l1l2seqs_75,[],1);

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
hist(gooddata.VYBmus(l1l2(scores>0.1)))
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)

setfig('mus');clf
hist(gooddata.VYBmus(l1l2))
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)

gooddata.loop1(l1l2(scores>0.1))
gooddata.loop2(l1l2(scores>0.1))

mean(gooddata.VYBmus(l1l2(scores>0.1)))
mean(gooddata.VYBmus(l1l2))

%% 
setfig('mu improvements');clf
endisA=l1l2seqs(:,end)==ones(size(l1l2'));
subplot(2,2,1)
title('End base is A');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisA))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisA))))
set(gca,'YScale','log')

subplot(2,2,2)
endisT=l1l2seqs(:,end)==ones(size(l1l2'))*2;

title('End base is T');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisT))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisT))))
set(gca,'YScale','log')

subplot(2,2,3)
endisC=l1l2seqs(:,end)==ones(size(l1l2'))*3;

title('End base is C');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisC))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisC))))
set(gca,'YScale','log')

subplot(2,2,4)
endisG=l1l2seqs(:,end)==ones(size(l1l2'))*4;

title('End base is G');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')
ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisG))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisG))))
set(gca,'YScale','log')
legend('All end bases', 'Fixed end base','location','best')


%% Starting base
%% 
setfig('mu improvements');clf
endisA=l1l2seqs(:,1)==ones(size(l1l2'));
subplot(2,2,1)
title('start base is A');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisA))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisA))))
set(gca,'YScale','log')

subplot(2,2,2)
endisT=l1l2seqs(:,1)==ones(size(l1l2'))*2;

title('start base is T');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisT))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisT))))
set(gca,'YScale','log')

subplot(2,2,3)
endisC=l1l2seqs(:,1)==ones(size(l1l2'))*3;

title('start base is C');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisC))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisC))))
set(gca,'YScale','log')

subplot(2,2,4)
endisG=l1l2seqs(:,1)==ones(size(l1l2'))*4;

title('start base is G');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')
ylabel('Count')

mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(endisG))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(endisG))))
set(gca,'YScale','log')
legend('All start bases', 'Fixed start base','location','best')

%% 
setfig('mu improvements others');clf
endisA=l1l2seqs(:,end)==ones(size(l1l2'));
startswithT=l1l2seqs(:,1)==ones(size(l1l2'))*2;
AT=endisA.*startswithT;

subplot(2,2,1)
title('End base is A, pairs with T');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

% mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(AT))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
mean(gooddata.VYBmus(l1l2(find(AT))))
set(gca,'YScale','log')

loop1endisA=l1l2seqs(:,6)==ones(size(l1l2'));
% startswithT=l1l2seqs(:,1)==ones(size(l1l2'))*2;
% AT=endisA.*startswithT;

subplot(2,2,2)
title('End base is A, pairs with T');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

% mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(loop1endisA))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
% mean(gooddata.VYBmus(l1l2(find(loop1endisA))))
set(gca,'YScale','log')




subplot(2,2,3)
tetraloop=(gooddata.loop1len==7).*( gooddata.loop2len==4) ;

title('tetraloop on loop1');
hold on
[n0,c0]=hist(gooddata.VYBmus,12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

% statistical test of significance of tetra and quintaloops
rp=randperm(length(gooddata.VYBmus),length(find(tetraloop)));
[h,p,ci,stats]=ttest(gooddata.VYBmus(find(tetraloop)),gooddata.VYBmus(rp))
mean(gooddata.VYBmus)
mean(gooddata.VYBmus(find(tetraloop)))

% mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(find(tetraloop)),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
% mean(gooddata.VYBmus(find(tetraloop)))
set(gca,'YScale','log')


subplot(2,2,4)
l1endG=l1l2seqs(:,6)==ones(size(l1l2'))*4;
l2startC=l1l2seqs(:,7)==ones(size(l1l2'))*3;

l1startT=l1l2seqs(:,1)==ones(size(l1l2'))*2;
l2endA=l1l2seqs(:,end)==ones(size(l1l2'))*1;

motif=(l1endG & l2startC) | (l1startT & l2endA);
motifATnGC=(l1endG & l2startC) & (l1startT & l2endA);
motifAT=(l1startT & l2endA);
title('End base is G');
hold on
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')
ylabel('Count')

% mean(gooddata.VYBmus(l1l2))

[n1,c1]=hist(gooddata.VYBmus(l1l2(find(motif))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
% mean(gooddata.VYBmus(l1l2(find(endisG))))
set(gca,'YScale','log')
legend('All end bases', 'Fixed end base','location','best')

mATi=find(motifAT);

rp=randperm(length(gooddata.VYBmus(l1l2(find(motifAT)))),length(find(motifATnGC)));
[h,p,ci,stats]=ttest(gooddata.VYBmus(l1l2(find(motifATnGC))),gooddata.VYBmus(l1l2(mATi(rp))))
mean(gooddata.VYBmus(l1l2(mATi)))
mean(gooddata.VYBmus(l1l2(find(motifATnGC))))
%% how does the number of A's affect mu
mATnGCi=find(motifATnGC);

setfig('As');clf
AAs=sum(l1l2seqs(mATnGCi,:)==ones(size(l1l2seqs(mATnGCi,:))),2);
plot(AAs,gooddata.VYBmus(l1l2(mATnGCi)),'.')

%% how can we get a string of A's... 
setfig('string of As');clf
l1twoAs=regexp(gooddata.loop1,'AAA+');
l1twoAsi=(~cellfun('isempty',l1twoAs));
shortloop1=gooddata.loop1len<6;

l2twoAs=regexp(gooddata.loop2,'AAA+');
l2twoAsi=(~cellfun('isempty',l2twoAs));
shortloop2=gooddata.loop2len<6;
motif=(l2twoAsi & shortloop2);

motif=(l1twoAsi & shortloop1);

hold on
[n0,c0]=hist(gooddata.VYBmus(shortloop2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')
ylabel('Count')

[n1,c1]=hist(gooddata.VYBmus(motif),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
% mean(gooddata.VYBmus(l1l2(find(endisG))))
set(gca,'YScale','log')


%% how can we continue looking for PK base pairs

motifs=zeros(length(l1l2seqs(:,5)),6);

leftbase=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];% A, T, T, G, C, G
rightbase=[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ];%T, A, G, T, G, C

left={'A', 'A', 'A', 'A', 'T', 'T','T','T', 'C','C','C','C','G','G','G','G'};
right={'A','T','C','G','A','T','C','G','A','T','C','G','A','T','C','G',};

% iterate 20 times
pvals=zeros(20,16);
for k=1:20
for i= 1:16
        l12ndlast=(l1l2seqs(:,5)==ones(size(l1l2'))*leftbase(i));
        l22ndstart=(l1l2seqs(:,8)==ones(size(l1l2'))*leftbase(i));
        motifs(:,i)= (l12ndlast & l22ndstart);
end

setfig('more pk');clf
for i=1:16
subplot(4,4,i)
[n0,c0]=hist(gooddata.VYBmus(l1l2),12);
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')
ylabel('Count')

hold on
[n1,c1]=hist(gooddata.VYBmus(l1l2(motifATnGC & motifs(:,i))),12);
area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
% mean(gooddata.VYBmus(l1l2(find(endisG))))
set(gca,'YScale','log')

sum(motifATnGC & motifs(:,i));
motifi=find(motifATnGC & motifs(:,i));

mean(gooddata.VYBmus(l1l2(mATnGCi)));
mean(gooddata.VYBmus(l1l2(motifATnGC & motifs(:,i))));

rp=randperm(length(gooddata.VYBmus(l1l2(mATnGCi))),length(motifi));
[h,p,ci,stats]=ttest(gooddata.VYBmus(l1l2(motifi)),gooddata.VYBmus(l1l2(mATnGCi(rp))));
text(-0.8,10^3,sprintf('p-value = %0.4f',p));
title(strcat(left(i),right(i)));
pvals(k,i)=p;
end
end



%% plot enrichment versus mu
setfig('enrichment mu correlation');clf
plot(gooddata.VYBmus,sum(gooddata.bincounts,2),'.')
xlabel('\mu')
ylabel('coverage')



%% END








