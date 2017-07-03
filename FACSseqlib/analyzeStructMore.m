addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
g=load('~/Documents/YFS/YFSI_gooddataWstruct.mat');
g=g.gooddata;

%% Plot distribution of fraction cleavers and mus
setfig('distrbution frac clv');clf
% plot(g.mus,g.fraccleavers,'.')
plot(g.mus,g.ensembstructfrac,'.')

%% Find most abundant structural motif
[c,ia,ic]=unique(g.ensembstruct);
% the following four lines is to get the occurrences of each unique
% string
y=sort(ic);
p = find([numel(y);diff(y);numel(y)]);
values = y(p(1:end-1));
instances = diff(p);

structmotif=c{instances==max(instances)};
o={};
for i=1:length(g.ensembstruct)
    o{end+1}=regexp(g.ensembstruct{i},structmotif,'tokens');
end

%%
structind=regexp(g.ensembstruct,structmotif);
setfig('structureless mus');clf
subplot(2,1,1)
hist(g.mus)
mean(g.mus)
subplot(2,1,2)
hist(g.mus(~cellfun('isempty',structind)))
mean(g.mus(~cellfun('isempty',structind)))
sum(~cellfun('isempty',structind))
% ensembstructfrac(end+1)=max(instances)/cnt;

%% What are the stem lengths?
longstem='(^\(+)(\.*)(\)+$)';
haslongstem1=~cellfun('isempty',regexp(g.loop1struct,longstem));
stem1length=[];
stem2length=[];
for i=1:length(g.loop1struct)
    o=regexp(g.loop1struct{i},longstem,'tokens');
    try
        stem1length(end+1)=length(strcat(o{1}{1},o{1}{3}))/2+6;
    catch
        stem1length(end+1)=6;
    end
end
fractionlongstem1=sum(haslongstem1)/length(haslongstem1); % 0.1694

for i=1:length(g.loop2struct)
    o=regexp(g.loop2struct{i},longstem,'tokens');
    try
        stem2length(end+1)=length(strcat(o{1}{1},o{1}{3}))/2+6;
    catch
        stem2length(end+1)=6;
    end
end

haslongstem2=~cellfun('isempty',regexp(g.loop2struct,longstem));
fractionlongstem2=sum(haslongstem2)/length(haslongstem2); % 0.0318

m1=zeros(max(stem1length),69);
mc1=m1;
m2=zeros(max(stem2length),69);
mc2=m2;
zerolengthloops=[];
loop=g.loop1len;

for k=1:length(loop)
    if loop(k)==0
        zerolengthloops(end+1)=k;
    else
        m1(stem1length(k),loop(k),end+1)=g.mus(k);
        mc1(stem1length(k),loop(k))=mc1(stem1length(k),loop(k))+1;
        m2(stem2length(k),loop(k),end+1)=g.mus(k);
        mc2(stem2length(k),loop(k))=mc2(stem2length(k),loop(k))+1;
    end
end

totaleachlength=zeros(1,69);
for k=1:69
    totaleachlength(k)=sum(loop==k);
end
mf1=mc1;
for i=1:max(stem1length)
    mf1(i,:)=mc1(i,:)./totaleachlength;
end

mf2=mc2;
for i=1:max(stem2length)
    mf2(i,:)=mc2(i,:)./totaleachlength;
end

%%
setfig('stem1 length');clf
imagesc(mf1)
load('MyColormaps','mycmap')
icmap=[[1.000 1.000 1.000];mycmap];
colormap(icmap)

set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 2 length')
ylabel('stem 1 length')
zlabel('frequency')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'frequency')

setfig('stem1 length mu dependence');clf
msum=sum(m1,3);
msum(msum==0)=NaN;
% mc(mc<10)=0;
morethan4=mc1>4;
mmean=msum./mf1./repmat(totaleachlength,max(stem1length),1).*morethan4;
imagesc(mmean)
load('MyColormaps','mycmap')
icmap=[[1.000 1.000 1.000];mycmap(end:-1:1,:)];
colormap(icmap)
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 2 length')
ylabel('stem 1 length')
zlabel('\mu')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'\mu')


setfig('stem1 length min mu dependence');clf
m1(m1==0)=NaN;
minm=min(m1,[],3);
% mc(mc<10)=0;
imagesc(minm)
load('MyColormaps','mycmap')
icmap=[[1.000 1.000 1.000];mycmap(end:-1:1,:)];
colormap(icmap)
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 2 length ')
ylabel('stem 1 length')
zlabel('\mu')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'min \mu')
grid on

%%
setfig('stem2 length');clf
imagesc(mf2)
load('MyColormaps','mycmap')
icmap=[[1.000 1.000 1.000];mycmap];
colormap(icmap)

set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 1 length')
ylabel('stem 2 length')
zlabel('frequency')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'frequency')

setfig('stem2 length mu dependence');clf
msum=sum(m2,3);
msum(msum==0)=NaN;
% mc(mc<10)=0;
morethan4=mc2>2;
mmean=msum./mf2./repmat(totaleachlength,max(stem2length),1).*morethan4;
imagesc(mmean)
load('MyColormaps','mycmap')
icmap=[[1.000 1.000 1.000];mycmap(end:-1:1,:)];
colormap(icmap)
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 1 length')
ylabel('stem 2 length')
zlabel('\mu')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'\mu')
%%

setfig('stem2 length min mu dependence');clf
m2(m2==0)=NaN;
minm=min(m2,[],3);
% mc(mc<10)=0;
imagesc(minm)
load('MyColormaps','mycmap')
icmap=[[1.000 1.000 1.000];mycmap(end:-1:1,:)];
colormap(icmap)
set(gca,'fontsize',14)
set(gca,'linewidth',1.5)
xlabel('loop 1 length')
ylabel('stem 2 length')
zlabel('\mu')
c=colorbar;
set(c,'fontsize',14)
set(c,'linewidth',1.5)
ylabel(c,'min \mu')
grid on

%% 
setfig('minimum vs count stem1');clf
plot(minm(:),mc1(:),'o','MarkerSize',5,'linewidth',2)
set(gca,'YScale','log')
set(gca,'linewidth',1.5)
set(gca,'fontsize',13)
xlabel('min \mu in category stem I')
ylabel('count')

setfig('minimum vs count stem2');clf
plot(minm(:),mc2(:),'o','MarkerSize',5,'linewidth',2)
set(gca,'YScale','log')
set(gca,'linewidth',1.5)
set(gca,'fontsize',13)
xlabel('min \mu in category')
ylabel('count')


%% 
setfig('mean vs count stem1');clf
plot(mmean(:),mc1(:),'o','MarkerSize',5,'linewidth',2)
set(gca,'YScale','log')
set(gca,'linewidth',1.5)
set(gca,'fontsize',13)
xlabel('mean \mu in stem II category')
ylabel('count')

setfig('mean vs count stem2');clf
plot(mmean(:),mc2(:),'o','MarkerSize',5,'linewidth',2)
set(gca,'YScale','log')
set(gca,'linewidth',1.5)
set(gca,'fontsize',13)
xlabel('mean \mu in stem II category')
ylabel('count')
%%
[stem1,loop1]=ind2sub(size(minm),find(minm==min(min(minm)))) 

%% are there sequences out there with similar structural motif as NL aptamer?
NLmotif='(^\.\.\.\.\.\.\.)(\(+.*\)+)(\.\.\.\.\.)(\(+.*\)+)(\.$)';
hasNL1=~cellfun('isempty',regexp(g.loop1struct,NLmotif)); % 1 total
hasNL2=~cellfun('isempty',regexp(g.loop2struct,NLmotif)); % 2 total

NLmus=g.mus([find(hasNL1) find(hasNL2)])
NLstruct=g.ensembstruct([find(hasNL1) find(hasNL2)])
NLseq=g.seqs{[find(hasNL1) find(hasNL2)]}
%%
nohairpin='^\(*\.+\)*$';
twoway='(^\(*\.+)(\(+.*\)+)(\.+\)*$)';
threeway='(^\(*\.+)(\(+.*\)+)(\.+)(\(+.*\)+)(\.+\)*$)';
fourway='(^\(*\.+)(\(+.*\)+)(\.+)(\(+.*\)+)(\.+)(\(+.*\)+)(\.+\)*$)';
fiveway='(^\(*\.+)(\(+.*\)+)(\.+)(\(+.*\)+)(\.+)(\(+.*\)+)(\.+)(\(+.*\)+)(\.+\)*$)'; %empty
loopstruct=g.loop2struct;

setfig('hist');clf
subplot(5,2,1)
hasnohairpin=~cellfun('isempty',regexp(loopstruct,nohairpin));
structnohairpin=g.ensembstruct(find(hasnohairpin));
seqnohairpin=g.seqs{find(hasnohairpin)};
musnohairpin=g.mus(find(hasnohairpin));
hist(musnohairpin)
axis([2 12 0 5000])

subplot(5,2,3)
has2way=~cellfun('isempty',regexp(loopstruct,twoway));
struct2way=g.ensembstruct(find(has2way));
seq2way=g.seqs{find(has2way)};
mus2way=g.mus(find(has2way));
hist(mus2way)
axis([2 12 0 5000])

subplot(5,2,5)
has3way=~cellfun('isempty',regexp(loopstruct,threeway));
struct3way=g.ensembstruct(find(has3way));
seq3way=g.seqs{find(has3way)};
mus3way=g.mus(find(has3way));
hist(mus3way)
axis([2 12 0 600])

subplot(5,2,7)
has4way=~cellfun('isempty',regexp(loopstruct,fourway));
struct4way=g.ensembstruct(find(has4way));
seq4way=g.seqs{find(has4way)};
mus4way=g.mus(find(has4way));
hist(mus4way)
axis([2 12 0 40])

subplot(5,2,9)
has5way=~cellfun('isempty',regexp(loopstruct,fiveway));
struct5way=g.ensembstruct(find(has5way));
seq5way=g.seqs{find(has5way)};
mus5way=g.mus(find(has4way));
hist(mus5way)
axis([2 12 0 35])

loopstruct=g.loop1struct;
subplot(5,2,2)
hasnohairpin=~cellfun('isempty',regexp(loopstruct,nohairpin));
structnohairpin=g.ensembstruct(find(hasnohairpin));
seqnohairpin=g.seqs{find(hasnohairpin)};
musnohairpin=g.mus(find(hasnohairpin));
hist(musnohairpin)
axis([2 12 0 6000])

subplot(5,2,4)
has2way=~cellfun('isempty',regexp(loopstruct,twoway));
struct2way=g.ensembstruct(find(has2way));
seq2way=g.seqs{find(has2way)};
mus2way=g.mus(find(has2way));
hist(mus2way)
axis([2 12 0 2000])

subplot(5,2,6)
has3way=~cellfun('isempty',regexp(loopstruct,threeway));
struct3way=g.ensembstruct(find(has3way));
seq3way=g.seqs{find(has3way)};
mus3way=g.mus(find(has3way));
hist(mus3way)
axis([2 12 0 300])

subplot(5,2,8)
has4way=~cellfun('isempty',regexp(loopstruct,fourway));
struct4way=g.ensembstruct(find(has4way));
seq4way=g.seqs{find(has4way)};
mus4way=g.mus(find(has4way));
hist(mus4way)
axis([2 12 0 20])




%% END



