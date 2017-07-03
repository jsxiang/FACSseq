addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS


%% load all data
yfsI=load('YFSI_parsed_countsonly19190.txt');
save('YFS19190.mat','yfsI')

%% plot distribution of barcodes
count_eachbarcode=sum(yfsI)/sum(sum(yfsI));
cellsorted=[1757865 2810 5817 31017 68874 257998 505385 1477277 2298559 1735153 698385 180782 1515006 2261 5166 24033 63119 221561 462551 1360497 2110373 2074741 483428 202671];
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
scaledcounts=yfsI.*repmat(scalefactor,length(yfsI(:,1)),1);
save('YFSI5252_scaled.mat','scaledcounts')
hold off
%% pick sequences with good statistics
rep1=1:12;
rep2=13:24;
tot1=sum(yfsI(:,rep1),2);
tot2=sum(yfsI(:,rep2),2);

rsquared=[];
num=[];
for n=10:100
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
plot(10:100,rsquared,'linewidth',1.5)
plot(10:100,num/length(goodstats),'linewidth',1.5)
xlabel('n')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
legend('r^2','fraction of seqs','location','best')
%%
n=44
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
rsq=mdl.Rsquared.Adjusted;

% plot(allmus(:,1),allmus(:,2),'.','MarkerSize',15)
dscatter(allmus(:,1),allmus(:,2),'BINS',[50 50])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('2-')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,n);
text(2,10,t)
%%
n=42
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
rsq=mdl.Rsquared.Adjusted;

% plot(allmus(:,1),allmus(:,2),'.','MarkerSize',15)
dscatter(allmus(:,1),allmus(:,2),'BINS',[50 50])
hold on
plot(x,0.92228*x+0.75506,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('2-')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,n);
text(2,10,t)

%% find histogram of distribution of coverage
setfig('coverage dist');clf
hist(tot1+tot2,40)
title('of 19190 sequences')
xlabel('coverage')
ylabel('no. of sequences')


%% combining both replicates and plot histogram of activity
goodscaledcountscombined=goodscaledcounts(:,rep1)+goodscaledcounts(:,rep2);
combinedmus=zeros(1,length(goodscaledcountscombined));

for i=1:length(goodscaledcounts(:,1))
    c1=goodscaledcounts(i,rep1);
    edges=1:12;
    [m,s]=normfit(edges,[],[],c1);
    combinedmus(i)=m;
end

setfig('hist mus');clf
hist(combinedmus,24)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('mu')
ylabel('count')

%% look at lengths
fid=fopen('YFSI_MoreThanOnce.txt');
seqs=textscan(fid,'%s');
fclose(fid);

yl=[];
backbonelength=length('GCTGTCACCGGATCCGGTCTGATGAGTCCGGACGAAACAGC');

for i=1:length(seqs{:})
    yl(end+1)=length(seqs{:}{i})-backbonelength;    
end


setfig('length');clf
hist(yl,60)
xlabel('lengths')
ylabel('count')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

%% look at lengths
fid=fopen('YFSI_MoreThanOnceShuf.txt');
seqs=textscan(fid,'%s');
fclose(fid);

yl=[];
backbonelength=length('GCTGTCACCGGATCCGGTCTGATGAGTCCGGACGAAACAGC');

for i=1:length(seqs{:})
    yl(end+1)=length(seqs{:}{i})-backbonelength;    
end


setfig('length');clf
hist(yl,60)
xlabel('lengths')
ylabel('count')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

%% look at lengths 
fid=fopen('YFSI_MoreThanOnceShuf19190.txt');
seqs=textscan(fid,'%s');
fclose(fid);

yl=[];
backbonelength=length('GCTGTCACCGGATCCGGTCTGATGAGTCCGGACGAAACAGC');

for i=1:length(seqs{:})
    yl(end+1)=length(seqs{:}{i})-backbonelength;    
end


setfig('length of 19190');clf
hist(yl,60)
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

%% plot each loop
fid=fopen('YFSI_19190loop1.txt');
loop1=textscan(fid,'%s');
fclose(fid);

fid=fopen('YFSI_19190loop2.txt');
loop2=textscan(fid,'%s');
fclose(fid);

l1len=[];
l2len=[];
for i=1:length(loop1{:})
    l1len(end+1)=length(loop1{:}{i});
end
for i=1:length(loop2{:})
    l2len(end+1)=length(loop2{:}{i});
end

i=1;
offset=[];
for k=1:length(yl)
    offset(end+1)=yl(k)-(l2len(k)+l1len(k));
end

%% Extract the loop sequences
% note, many of the sequences are chimeras, we should throw them out
o=regexp(s,'(GCTGTCACCGGA)([A|C|T|G]*)(TCCGGTCTGATGAGTCC)','tokens')








%% END 

