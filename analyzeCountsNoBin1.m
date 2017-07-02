addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/Standard/

addpath /Volumes/smolke-lab$/Joy/YFSI/
clear
numseqs=15
% numseqs=1
% from validation results, m= 0.2156  , c=-1.6998, for x=NGSmu
validm=0.2156;
validc=-1.6998;

%% load all data
yfsI=load('YFSI_parsed_uniqcountsonly.txt');
% save('YFSI.mat','yfsI');

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
save('YFSI_scaled.mat','scaledcounts')
hold off



%% pick sequences with good statistics
% rep1=2:12;
% rep2=14:24;
% tot1=sum(yfsI(:,rep1),2);
% tot2=sum(yfsI(:,rep2),2);
% 
% rsquared=[];
% num=[];
% for n=10:2:250
%     goodstats=tot1>n & tot2>n;
% 
%     goodscaledcounts=scaledcounts(find(goodstats),:);
% 
%     %fit data to normdist
% 
%     rep1mus=zeros(1,length(goodscaledcounts(:,1)));
%     rep2mus=rep1mus;
%     for i=1:length(goodscaledcounts(:,1))
%         c1=goodscaledcounts(i,rep1);
%         edges=2:12;
%         fedges=-1.6998+0.2156*edges;
%         [m,s]=normfit(edges,[],[],c1);
%         rep1mus(i)=m;
% 
%         c2=goodscaledcounts(i,rep2);
%         [m,s]=normfit(edges,[],[],c2);
%         rep2mus(i)=m;
%     end
% 
%     % plot
% 
%     setfig('compare replicates');clf
%     x=linspace(min(edges),max(edges));
%     y=x;
% 
%     allmus=[rep1mus' rep2mus'];
% 
%     mdl = fitlm(allmus(:,1),allmus(:,2));
%     rsquared(end+1)=mdl.Rsquared.Adjusted;
%     num(end+1)=sum(goodstats);
% end
% 
% setfig('n versus rsquared');clf
% hold on
% plot(10:2:250,rsquared,'linewidth',1.5)
% plot(10:2:250,num/length(goodstats),'linewidth',1.5)
% xlabel('n')
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% legend('r^2','fraction of seqs','location','best')

%%
rep1=2:12;
rep2=14:24;
tot1=sum(yfsI(:,rep1),2);
tot2=sum(yfsI(:,rep2),2);

goodstats=tot1>numseqs & tot2>numseqs;
goodscaledcounts=scaledcounts(find(goodstats),:);
origcounts=yfsI(find(goodstats),:);

%fit data to normdist
rep1mus=zeros(1,length(goodscaledcounts(:,1)));
rep2mus=rep1mus;
for i=1:length(goodscaledcounts(:,1))
    c1=goodscaledcounts(i,rep1);
    edges=2:12;
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
xlabel('\mu rep 1')
ylabel('\mu rep 2')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn > %0.0f',rsq,numseqs);
text(-0.8,5,t)

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
    edges=2:12;
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
fid=fopen('YFSI_parsed_uniqcountsseqs.txt');
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
gooddata.VYBmus=validm*gooddata.mus+validc;
gooddata.VYBmusall=validm*allmus+validc;

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

save('YFSI_gooddata_noBin1.mat','gooddata');



%% plot with VYBmus
setfig('compare replicates vybmus');clf
x=linspace(min(min(gooddata.VYBmusall)),max(max(gooddata.VYBmusall)));
y=x;

mdl = fitlm(gooddata.VYBmusall(:,1),gooddata.VYBmusall(:,2));
rsq=mdl.Rsquared.Adjusted;

% plot(allmus(:,1),allmus(:,2),'.','MarkerSize',15)
dscatter(gooddata.VYBmusall(:,1),gooddata.VYBmusall(:,2),'BINS',[50 50])
hold on
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
hold off
xlabel('\mu rep 1')
ylabel('\mu rep 2')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn > %0.0f',rsq,numseqs);
text(-0.8,0.5,t)
axis([-1.2 1 -1.2 1])
%% Pick 12 spaced out activity sequences to test
setfig('pick');clf
dscatter(allmus(:,1),allmus(:,2),'BINS',[50 50])
goodagreement=find(abs(mdl.Coefficients.Estimate(2)*rep1mus+mdl.Coefficients.Estimate(1)-rep2mus)<0.4);
hold on
t=sprintf('R^2 = %0.2f\nn > %0.0f',rsq,numseqs);
text(5,10,t)
plot(rep1mus(goodagreement),rep2mus(goodagreement),'r.','MarkerSize',10)
% plot(combinedmus(goodagreement),combinedmus(goodagreement),'oy','MarkerSize',20)
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlabel('rep 1 \mu')
ylabel('rep 2 \mu')

murange=max(combinedmus(goodagreement))-min(combinedmus(goodagreement));
mucand=linspace(0.01*murange+min(combinedmus(goodagreement)),0.95*murange+min(combinedmus(goodagreement)),12);
mucand=linspace(0.01*murange+min(combinedmus(goodagreement)),0.10*murange+min(combinedmus(goodagreement)),12);
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


%% 
f=fopen('goodseqs.txt','w');
fprintf(f,'%s\n',gooddata.seqs{:});
fclose(f);

%%
gooddata.standmus=Standard(gooddata.VYBmus,2,mean(gooddata.VYBmus),std(gooddata.VYBmus));

g=fopen('goodmus.txt','w');
fprintf(g,'%0.4f\n',gooddata.standmus);
fclose(g);


%% go through all sequences and convert to matrix of 1 2 3 4
% note: full length sequences are considered
seq=cell(length(gooddata.seqs),1);
% loop2seq=loop1seq;
for i=1:length(gooddata.seqs)
    s=gooddata.loop1{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop1seq{i}=s-'0';
%     
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
%%
totlen=113;
l1l2seqs=zeros(length(gooddata.seqs),totlen);

for i=1:length(gooddata.seqs)
    l1l2seqs(i,13:(12+length(gooddata.loop1seqnum{1,i})))=gooddata.loop1seqnum{1,i};
    l1l2seqs(i,(12+length(gooddata.loop1seqnum{1,i})+17+1):(12+length(gooddata.loop1seqnum{1,i})+17+length(gooddata.loop2seqnum{1,i})))=gooddata.loop2seqnum{1,i};
end

%% convert to bit 
bitmat1=zeros(length(gooddata.seqs), totlen,4);

for i=1:length(gooddata.seqs)
    for j=1:4
        for k=1:totlen
            bitmat1(i,k,j)=l1l2seqs(i,k)==j;
        end
    end
end

%% 
setfig('enrichment vs. mu');clf
plot(gooddata.VYBmus,sum(gooddata.bincounts,2),'.')
hold on
plot(gooddata.VYBmus,sum(gooddata.origcounts,2),'.')
set(gca,'YScale','log')
xlabel('\mu')
ylabel('counts')
prctile(gooddata.VYBmus,1)
%% 
f=fopen('goodseqs_onehot.txt','w');
fprintf(f,'%i %i\n',bitmat1(1:10,1:2));
fclose(f);

%% Train, valid, test
ri=randperm(length(bitmat1),length(bitmat1));
bitmatrand=bitmat1(ri,:,:);
switches=gooddata.VYBmus<prctile(gooddata.VYBmus,15);

switchesrand=switches(ri)';
tr=struct;
tr.trainxdata=bitmatrand(1:round(0.8*length(bitmatrand(:,1,1))),:,:);
tr.traindata=switchesrand(1:round(0.8*length(bitmatrand(:,1,1))));
save('trainALL.mat','tr');

v=struct;
v.validxdata=bitmatrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))),:,:);
v.validdata=switchesrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))));
save('validALL.mat','v');

tt=struct;
tt.testxdata=bitmatrand((round(0.9*length(bitmatrand(:,1,1)))+1):end,:,:);
tt.testdata=switchesrand((round(0.9*length(bitmatrand(:,1,1)))+1):end);
save('testALL.mat','tt');












%% END 

