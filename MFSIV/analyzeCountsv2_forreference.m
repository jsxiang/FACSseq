addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
clear
numseqs=15

%% load all data
yfsII=load('YFSII_parsed_uniqcountsonly.txt');
% save('YFSI.mat','yfsI');

%% plot distribution of barcodes
count_eachbarcode=sum(yfsII)/sum(sum(yfsII));
cellsorted=[87296 87296 698368 698368 43008 43008 519460 689166 2748613 2878572 2131946 2008886 1551192 1498790 507997 403930 176413 173248 73001 57827 25872 24072 10617 8401 3062 3518 2425 1778 1186 1598];
cellsorted=cellsorted/sum(cellsorted);
setfig('barcode distribution');clf
plot(cellsorted,count_eachbarcode,'.','MarkerSize',15)
xlabel('cells sorted')
ylabel('NGS count')

hold on
xval=linspace(min(cellsorted),max(cellsorted));
yval=xval;
plot(xval,yval)
ylow=xval-0.02;
yhigh=xval+0.02;
plot(xval,ylow,xval,yhigh)

barcodestoohigh=find(count_eachbarcode>(cellsorted+0.02));
barcodestoolow=find(count_eachbarcode<(cellsorted-0.02));
scalefactor=cellsorted./count_eachbarcode;
scaledcounts=yfsII.*repmat(scalefactor,length(yfsII(:,1)),1);
save('YFSII_scaled.mat','scaledcounts')
hold off
%% pick sequences with good statistics
rep1=1:12;
rep2=13:24;
tot1=sum(yfsII(:,rep1),2);
tot2=sum(yfsII(:,rep2),2);

rsquared=[];
num=[];
for n=10:2:250
    goodstats=tot1>n & tot2>n;

    goodscaledcounts=scaledcounts(find(goodstats),:);

    %fit data to normdist

    rep1mus=zeros(1,length(goodscaledcounts(:,1)));
    rep2mus=rep1mus;
    for i=1:length(goodscaledcounts(:,1))
        c1=goodscaledcounts(i,rep1);
        edges=1:12;
        [m,s]=normfit(edges,[],[],c1);
        rep1mus(i)=m;

        c2=goodscaledcounts(i,rep2);
        [m,s]=normfit(edges,[],[],c2);
        rep2mus(i)=m;
    end

    % plot

    setfig('compare replicates');clf
    x=linspace(min(edges),max(edges));
    y=x;

    allmus=[rep1mus' rep2mus'];

    mdl = fitlm(allmus(:,1),allmus(:,2));
    rsquared(end+1)=mdl.Rsquared.Adjusted;
    num(end+1)=sum(goodstats);
end

setfig('n versus rsquared');clf
hold on
plot(10:2:250,rsquared,'linewidth',1.5)
plot(10:2:250,num/length(goodstats),'linewidth',1.5)
xlabel('n')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
legend('r^2','fraction of seqs','location','best')

%%
rep1=1:12;
rep2=13:24;
tot1=sum(yfsII(:,rep1),2);
tot2=sum(yfsII(:,rep2),2);

goodstats=tot1>numseqs & tot2>numseqs;
goodscaledcounts=scaledcounts(find(goodstats),:);
origcounts=yfsII(find(goodstats),:);

%fit data to normdist
rep1mus=zeros(1,length(goodscaledcounts(:,1)));
rep2mus=rep1mus;
for i=1:length(goodscaledcounts(:,1))
    c1=goodscaledcounts(i,rep1);
    edges=1:12;
    [m,s]=normfit(edges,[],[],c1);
    rep1mus(i)=m;

    c2=goodscaledcounts(i,rep2);
    [m,s]=normfit(edges,[],[],c2);
    rep2mus(i)=m;
end

% plot
setfig('compare replicates');clf
x=linspace(min(edges),max(edges));
y=x;

allmus=[rep1mus' rep2mus'];

mdl = fitlm(allmus(:,1),allmus(:,2));
rsq=mdl.Rsquared.Adjusted;

% plot(allmus(:,1),allmus(:,2),'.','MarkerSize',15)
dscatter(allmus(:,1),allmus(:,2),'BINS',[50 50])
hold on
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('2-')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,numseqs);
text(2,10,t)

%% find histogram of distribution of coverage
setfig('coverage dist');clf
hist(tot1+tot2,500)
title('of all sequences')
xlabel('coverage')
ylabel('no. of sequences')


%% combining both replicates and plot histogram of activity
goodscaledcountscombined=goodscaledcounts(:,rep1)+goodscaledcounts(:,rep2);
combinedmus=zeros(1,length(goodscaledcountscombined));
combinedsigmas=combinedmus;
for i=1:length(goodscaledcounts(:,1))
    c1=goodscaledcounts(i,rep1)+goodscaledcounts(i,rep2);
    edges=1:12;
    [m,s]=normfit(edges,[],[],c1);
    combinedmus(i)=m;
    combinedsigmas(i)=s;
end

setfig('hist mus');clf
hist(combinedmus,12)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('mu')
ylabel('count')

%% extract sequences
fid=fopen('YFSII_parsed_uniqcountsseqs.txt');
seqs=textscan(fid,'%s');
fclose(fid);

yl=[];
backbonelength=length('GCTGTCACCGGATCCGGTCTGATGAGTCCGGACGAAACAGC');

for i=1:length(seqs{:})
    yl(end+1)=length(seqs{:}{i})-backbonelength;    
end

setfig('length');clf
hist(yl,max(yl)-min(yl))
xlabel('lengths')
ylabel('count')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

%% plot lengths with activity
setfig('lengths vs mu');clf
dscatter(yl(goodstats)',combinedmus')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('lengths')
ylabel('\mu')


%% Extract the loop sequences
% note, many of the sequences are chimeras, we should throw them out
goodind=find(goodstats);
goodseqs={};
for i=1:sum(goodstats) 
    goodseqs{end+1}=seqs{1}{goodind(i)};
end

loop1={};
loop2={};
loop1len=[];
loop2len=[];
for i=1:length(goodseqs)
    o=regexp(goodseqs{i},'(GCTGTCACCGGA)([A|C|T|G]*)(TCCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)','tokens');
    loop1{end+1}=o{1}{2};
    loop2{end+1}=o{1}{4};
    loop1len(end+1)=length(o{1}{2});
    loop2len(end+1)=length(o{1}{4});
end

gooddata=struct;
gooddata.seqs=goodseqs;
gooddata.bincounts=goodscaledcounts;
gooddata.origcounts=origcounts;
gooddata.mus=combinedmus;
gooddata.sigmas=combinedsigmas;
gooddata.loop1=loop1;
gooddata.loop2=loop2;
gooddata.loop1len=loop1len;
gooddata.loop2len=loop2len;
gindex={};
gi=find(goodstats);
for i=1:length(gi)
    ind=sprintf('YI%0.0f',gi(i));
    gindex{end+1}=ind;
end
gooddata.index=gindex;
gooddata.indexnum=gi;
    


setfig('loop1 length versus loop2 length');clf
dscatter(gooddata.loop1len',gooddata.loop2len','bins',[69 69],'marker','s')
xlabel('loop1 length')
ylabel('loop2 length')
set(gca,'FontSize',14)
set(gca,'linewidth',1.5)

save('YFSI_gooddata.mat','gooddata');



%% Pick 12 spaced out activity sequences to test
setfig('pick');clf
dscatter(allmus(:,1),allmus(:,2),'BINS',[50 50])
goodagreement=find(abs(mdl.Coefficients.Estimate(2)*rep1mus+mdl.Coefficients.Estimate(1)-rep2mus)<0.4);
hold on
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,numseqs);
text(2,10,t)
plot(rep1mus(goodagreement),rep2mus(goodagreement),'r.','MarkerSize',10)
% plot(combinedmus(goodagreement),combinedmus(goodagreement),'oy','MarkerSize',20)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('rep 1 \mu')
ylabel('rep 2 \mu')

murange=max(combinedmus(goodagreement))-min(combinedmus(goodagreement));
mucand=linspace(0.025*murange,0.95*murange,12);
mucandind=[];
for i=1:length(mucand)
    candind=find(abs(mucand(i)-combinedmus(goodagreement))==min(abs(mucand(i)-combinedmus(goodagreement))));    
    mucandind(end+1)=candind(1);
        
end
plot(rep1mus(goodagreement(mucandind)),rep2mus(goodagreement(mucandind)),'ko','MarkerSize',10,'linewidth',2)
setfig('where are the mus');clf
bar([rep1mus(goodagreement(mucandind));rep2mus(goodagreement(mucandind))]')
seqcand=gooddata.seqs(goodagreement(mucandind))
%% check that these candidate sequences are validate sequences
clc
for i=1:length(mucandind)
    s=gooddata.seqs{goodagreement(mucandind(i))};
    o=regexp(s,'(GCTGTCACCGGA)([A|C|T|G]*)(TCCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)','tokens');
    fprintf('%s\t%s\n',o{1}{[2,4]})
end
for i=1:length(mucandind)
    s=gooddata.seqs{goodagreement(mucandind(i))};
    fprintf('%s\n',s)
end
for i=1:length(mucandind)
    s=gooddata.index{goodagreement(mucandind(i))};
    fprintf('%s\n',s)
end
%% check that the distributions of the bin counts look reasonable as well
setfig('check bin counts');clf
for i=1:length(mucandind)
    subplot(12,2,2*i-1)
    bc=gooddata.bincounts(goodagreement(mucandind(i)),rep1);
    bar(bc)
    title('rep 1')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',10)
    xlabel('bin')
    ylabel('count')
    
    subplot(12,2,2*i)
    bc=gooddata.bincounts(goodagreement(mucandind(i)),rep2);
    bar(bc)    
    title('rep 2')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',10) 
    xlabel('bin')
    ylabel('count')
end
%% END 

