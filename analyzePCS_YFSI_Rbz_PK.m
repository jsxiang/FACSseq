addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/FACSseq/
% addpath ~/Volumes/smolke-lab$/Joy/YFSI/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
g=load('YFSI_gooddata_noBin1.mat');
gooddata=g.gooddata;


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
% l1len=21;
l1len=6;

endslen=0;
% totallooplength=l1len+l2len;
% l1l2=find(gooddata.loop2len==l2len);
l1l2=find((gooddata.loop1len==l1len)'.* (gooddata.VYBmus(:,1)<prctile(gooddata.VYBmus(:,1),25)));

l1l2seqs=zeros(length(l1l2),l1len);
for i=1:length(l1l2)
    l1l2seqs(i,1:l1len)=gooddata.loop1seqnum{l1l2(i)};
%     l1l2seqs(i,1:(l2len))=gooddata.loop2seqnum{l1l2(i)};
end


%% f
setfig('mu dist');clf
[n0,c0]=hist(gooddata.VYBmus(:,1),20)
[n1,c1]=hist(gooddata.VYBmus(:,2),20)
hold on
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
set(gca,'YScale','log')
hold off
legend('no FA','with 1mM FA','location','best')

%% Analyses

muth=gooddata.VYBmus<prctile(gooddata.VYBmus,100);
l2len=[15:25 30 40 50 60];
l1len=[5 6];
endslen=[3 3]
[l1l2seqs,l1l2]=findSeqs(gooddata,l1len,l2len,endslen,muth);

% what is the distribution of mus
findParameters(gooddata.VYBmus,l1l2seqs);

% find mutual sequence contribution - pairwise mu
pairwise=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,15,1);

% find entropy
ent=findEnt(l1l2seqs);

% find mutual information
[mI,pJoint,pointMI]=findMutualInformation(l1l2seqs);

%% perform consensus analysis with pairwise mu
pair=pairwise;
pair(find(pair==min(min(pair))))=nan;
[row,col]=find(pair==min(min(pair)));

base1=rem(row(1)-1,4)+1
pos1=ceil(row(1)/4)

loop2re=regexp(gooddata.loop2,'^G');
loop2consensus1=(~cellfun('isempty',loop2re));

base2=rem(row(2)-1,4)+1
pos2=ceil(row(2)/4)
loop1re=regexp(gooddata.loop1,'^G');
loop1consensus1=(~cellfun('isempty',loop1re));

%%
muth=loop1consensus1.*loop2consensus1;

[l1l2seqs,l1l2]=findSeqs(gooddata,l1len,l2len,endslen,muth);

pairwise=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,15,1);
pair=pairwise;
pair(find(pair==min(min(pair))))=nan;
[row,col]=find(pair==min(min(pair)));
base1=rem(row(1)-1,4)+1
pos1=ceil(row(1)/4)
base2=rem(row(2)-1,4)+1
pos2=ceil(row(2)/4)

loop2re=regexp(gooddata.loop1,'G$');
loop2consensus2=(~cellfun('isempty',loop2re));
loop1re=regexp(gooddata.loop1,'A[A|C|T|G][A|C|T|G]$');
loop1consensus2=(~cellfun('isempty',loop1re));
%%
muth=loop1consensus1.*loop2consensus1.*loop2consensus2.*loop1consensus2;
[l1l2seqs,l1l2]=findSeqs(gooddata,l1len,l2len,endslen,muth);

pairwise=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,15,1);
pair=pairwise;
pair(find(pair==min(min(pair))))=nan;
[row,col]=find(pair==min(min(pair)));
base1=rem(row(1)-1,4)+1
pos1=ceil(row(1)/4)
base2=rem(row(2)-1,4)+1
pos2=ceil(row(2)/4)

loop2re=regexp(gooddata.loop2,'^[A|T|C|G]G');
loop2consensus3=(~cellfun('isempty',loop2re));
loop1re=regexp(gooddata.loop1,'A[A|C|T|G]$');
loop1consensus3=(~cellfun('isempty',loop1re));

%%
muth=muth.*loop1consensus3.*loop2consensus3;
[l1l2seqs,l1l2]=findSeqs(gooddata,l1len,l2len,endslen,muth);

pairwise=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,15,1);
pair=pairwise;
pair(find(pair==min(min(pair))))=nan;
[row,col]=find(pair==min(min(pair)));
base1=rem(row(1)-1,4)+1
pos1=ceil(row(1)/4)
base2=rem(row(2)-1,4)+1
pos2=ceil(row(2)/4)

%%

muth=gooddata.VYBmus<prctile(gooddata.VYBmus,25);
[l1l2seqs_25,l1l2_25]=findSeqs(gooddata,l1len,l2len,endslen,muth);

muth=gooddata.VYBmus>prctile(gooddata.VYBmus,75);
[l1l2seqs_75,l1l2_75]=findSeqs(gooddata,l1len,l2len,endslen,muth);

% find mutual information
[mI_25,pJoint_25,allI_25]=findMutualInformation(l1l2seqs_25);

[mI_75,pJoint_75,allI_75]=findMutualInformation(l1l2seqs_75);



%% look for average length dependence of mu
m=zeros(69,69);
mc=m;
zerolengthloops=[];
for k=1:length(gooddata.seqs)
    if gooddata.loop1len(k)==0||gooddata.loop2len(k)==0
        zerolengthloops(end+1)=k;
    else
        m(gooddata.loop1len(k),gooddata.loop2len(k),end+1)=gooddata.mus(k);
        mc(gooddata.loop1len(k),gooddata.loop2len(k))=mc(gooddata.loop1len(k),gooddata.loop2len(k))+1;
    end
end
msum=sum(m,3);
msum(msum==0)=NaN;
mc(mc<10)=0;
mmean=msum./mc;
setfig('loop length dependence');clf
imagesc(mmean)
load('MyColormaps','mycmap')
icmap=mycmap(end:-1:1,:);
icmap=[[1.000 1.000 1.000];icmap];
colormap(icmap)
axis([0 70 0 70])

set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 2 length')
ylabel('loop 1 length')
zlabel('\mu')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'\mu')

setfig('number per length');clf
surf(mc)
mycmap=[[1.000 1.000 1.000];mycmap];
colormap(mycmap)
% c=colorbar;
% set(c,'fontsize',14)
% set(c,'linewidth',1.5)
% ylabel(c,'count')
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 2 length')
ylabel('loop 1 length')
zlabel('no. of seqs')
axis([0 70 0 70 0 max(max(max(mc)))])


%% look for sequence motifs in loop1

loop1motif='^T';
seqRE1=regexp(gooddata.loop1,loop1motif);
loop1m=(~cellfun('isempty',seqRE1));

loop2motif='A$';
seqRE2=regexp(gooddata.loop2,loop2motif);
loop2m=(~cellfun('isempty',seqRE2));

pk = find((loop1m.*loop2m));


% shortloop1=gooddata.loop1len<10;

setfig('pseudoknots');clf

[n0,c0]=hist(gooddata.VYBmus,20);
[n1,c1]=hist(gooddata.VYBmus(pk),20);
hold on
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
set(gca,'YScale','log')
hold off
legend('all short loop I',strcat('with pseudoknot',': ',loop1motif),'location','best')


map_i=containers.Map(gooddata.loop1,1:length(gooddata.loop1));

cind=[];
for i=1:length(loop1m)
    try
    cind(end+1)=map_i(gooddata.loop1{loop1m(i)});
    end
   
end
%%
x=linspace(min(gooddata.VYBmus(:,1)),max(gooddata.VYBmus(:,2)));

setfig('hits pk on loop1');clf
hold on
dscatter(gooddata.VYBmusall,gooddata.VYBmusall,'BINS',[50 50])
plot(x,x,'k:','linewidth',2.5)
plot(gooddata.VYBmusall(cind,1),gooddata.VYBmusall(cind,2),'ro','MarkerSize',10,'linewidth',2)
xlabel('rep 1')
ylabel('rep 2')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
title('\mu')
hold off

%% look for sequence motifs in loop2

loop2motif='G$';
seqRE2=regexp(gooddata.loop2,loop2motif);
loop2m=find(~cellfun('isempty',seqRE2));
shortloop2=gooddata.loop2len<10;

setfig('pseudoknots on loop2');clf

[n0,c0]=hist(gooddata.cmus(shortloop2,1),20);
[n1,c1]=hist(gooddata.cmus(loop2m,1),20);
hold on
area(c0,(n0),'facecolor',[1.0 0.8 0],'edgecolor',[1.0 0.6 0],'facealpha',0.2,'linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu')

ylabel('Count')

area(c1,(n1),'facecolor',[1.0 0.3 .3],'edgecolor',[1.0 0.3 .3],'facealpha',0.5,'linewidth',2)
set(gca,'YScale','log')
hold off
legend('all short loop II',strcat('with pseudoknot',': ',loop2motif),'location','best')


map_i=containers.Map(gooddata.loop2,1:length(gooddata.loop2));

cind=[];
for i=1:length(loop2m)
    try
    cind(end+1)=map_i(gooddata.loop2{loop2m(i)});
    end
   
end

x=linspace(min(gooddata.cmus(:,1)),max(gooddata.cmus(:,2)));

setfig('hits pk on loop2');clf
hold on
dscatter(gooddata.cmus(shortloop2,1),gooddata.cmus(shortloop2,2),'BINS',[50 50])
plot(x,x,'k:','linewidth',2.5)
plot(gooddata.cmus(cind,1),gooddata.cmus(cind,2),'ro','MarkerSize',10,'linewidth',2)
xlabel('0 mM Norlaudanosoline')
ylabel('1 mM Norlaudanosoline')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
title('\mu')
hold off
%%
rep1minus=29:-2:7;
rep1plus=30:-2:8;
setfig('check bin counts again');clf
for i=1:length(cind)
    subplot(length(cind),4,4*i-3)
    bc=gooddata.bincounts1(cind(i),rep1minus);
    bar(bc)
    title('rep 1-')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',10)
    xlabel('bin')
    ylabel('count')
    
    subplot(length(cind),4,4*i-2)
    bc=gooddata.bincounts1(cind(i),rep1plus);
    bar(bc)    
    title('rep 1+')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',10) 
    xlabel('bin')
    ylabel('count')
    
    subplot(length(cind),4,4*i-1)
    bc=gooddata.bincounts2(cind(i),rep1minus);
    bar(bc)
    title('rep 2-')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',10)
    xlabel('bin')
    ylabel('count')
    
    subplot(length(cind),4,4*i)
    bc=gooddata.bincounts2(cind(i),rep1plus);
    bar(bc)    
    title('rep 2+')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',10) 
    xlabel('bin')
    ylabel('count')
end





%% END








