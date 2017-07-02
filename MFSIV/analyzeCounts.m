addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/
clear
%%
numseqs=10;
findONswitches=0;
findOFFswitches=0;
findnonswitches=1;

% linearc=-3.2214;
% linearm=0.3889;

%% load  data
mfsIV=load('MFSIV_parsed_uniqcountsonly.txt');
% %% load all data
mfsIV=mfsIV(:,28:54);
%% plot distribution of barcodes
% count_eachbarcode=sum(mfsIV)/sum(sum(mfsIV));
count_eachbarcode=sum(mfsIV);
cellsorted=[193791      198832      196483      197448      350994      350954      381104      401279      243437      238285      198975      206054       88712       87346 75943       74949       31840       30570       32826       30339        7667        7092        9278        8769        6656        4608        5000];
% cellsorted=cellsorted/sum(cellsorted);
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
%%
barcodestoohigh=find(count_eachbarcode>(cellsorted+0.02));
barcodestoolow=find(count_eachbarcode<(cellsorted-0.02));
scalefactor=cellsorted./count_eachbarcode;
scaledcounts=mfsIV.*repmat(scalefactor,length(mfsIV(:,1)),1); % actually may not want to scale so early
% scaledcounts=mfsIV;
save('MFSIV_scaled.mat','scaledcounts')
hold off
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
title('MFSIV')

%% set problematic ones to zero
% problematic=[6 8 9 19 21 23];
% 
% scaledcounts(:,problematic)=0;
%%
clear yfsIII
clear yfsII


%% pick sequences with good statistics
rep1minus=21:-4:1;
rep1plus=23:-4:3;
tot1minus=sum(mfsIV(:,rep1minus),2);
tot1plus=sum(mfsIV(:,rep1plus),2);

rep2minus=22:-4:2;
rep2plus=24:-4:4;
tot2minus=sum(mfsIV(:,rep2minus),2);
tot2plus=sum(mfsIV(:,rep2plus),2);


% rsquared=[];
% num=[];
% for n=2:2:250
%     goodstats=tot1minus>n & tot1plus>n & tot2minus>n & tot2plus>n;
% 
%     goodscaledcounts=scaledcounts(find(goodstats),:);
% 
%     %fit data to normdist
% 
%     rep1minusmus=zeros(1,length(goodscaledcounts(:,1)));
%     rep1plusmus=rep1minusmus;
%     rep2minusmus=rep1minusmus;
%     rep2plusmus=rep1minusmus;    
%     for i=1:length(goodscaledcounts(:,1))
%         c1=goodscaledcounts(i,rep1minus);
%         edges=1:6;
%         [m,s]=normfit(edges,[],[],c1);
%         rep1minusmus(i)=m;
% 
%         c2=goodscaledcounts(i,rep1plus);
%         [m,s]=normfit(edges,[],[],c2);
%         rep1plusmus(i)=m;
%         
%         c1=goodscaledcounts(i,rep2minus);
%         [m,s]=normfit(edges,[],[],c1);
%         rep2minusmus(i)=m;
% 
%         c2=goodscaledcounts(i,rep2plus);
%         [m,s]=normfit(edges,[],[],c2);
%         rep2plusmus(i)=m;        
%     end
% 
% 
% 
%     allmus=[rep1minusmus' rep2minusmus' rep1plusmus' rep2plusmus'];
% 
%     mdl = fitlm(allmus(:,1),allmus(:,2));
%     rsquared(end+1)=mdl.Rsquared.Adjusted;
%     num(end+1)=sum(goodstats);
% end
% 
% x=linspace(min(edges),max(edges));
% y=x;
% 
% setfig('n versus rsquared');clf
% hold on
% plot(2:2:250,rsquared,'linewidth',1.5)
% plot(2:2:250,num/length(goodstats),'linewidth',1.5)
% xlabel('n')
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% legend('r^2','fraction of seqs','location','best')
% title('mfsIV')

%% Choose a numseq, plot data 

numseqs=25;
goodstats=tot1minus>numseqs & tot1plus>numseqs & tot2minus>numseqs & tot2plus>numseqs;
goodscaledcounts1=scaledcounts(find(goodstats),:);
origcounts1=mfsIV(find(goodstats),:);

%fit data to normdist
rep1minusmus=zeros(1,length(goodscaledcounts1(:,1)));
rep1plusmus=rep1minusmus;
rep1minussigs=rep1minusmus;
rep1plussigs=rep1minusmus;

rep2minusmus=zeros(1,length(goodscaledcounts1(:,1)));
rep2plusmus=rep1minusmus;
rep2minussigs=rep1minusmus;
rep2plussigs=rep1minusmus;

rep1mci=rep1minusmus;
rep2mci=rep1minusmus;
rep1sci=rep1minusmus;
rep2sci=rep1minusmus;
for i=1:length(goodscaledcounts1(:,1))
    c1=goodscaledcounts1(i,rep1minus);
    edges=1:6;
    [m,s,mci,sci]=normfit(edges,[],[],c1);
    rep1minusmus(i)=m;
    rep1minussigs(i)=s;
%     rep1mci(i)=mci;
%     rep1sci(i)=sci;
    
    c2=goodscaledcounts1(i,rep1plus);
    [m,s]=normfit(edges,[],[],c2);
    rep1plusmus(i)=m;
    rep1plussigs(i)=s;
%     rep2mci(i)=mci;
%     rep2sci(i)=sci;

    c1=goodscaledcounts1(i,rep2minus);
    [m,s,mci,sci]=normfit(edges,[],[],c1);
    rep2minusmus(i)=m;
    rep2minussigs(i)=s;
%     rep1mci(i)=mci;
%     rep1sci(i)=sci;
    
    c2=goodscaledcounts1(i,rep2plus);
    [m,s]=normfit(edges,[],[],c2);
    rep2plusmus(i)=m;
    rep2plussigs(i)=s;
%     rep2mci(i)=mci;
%     rep2sci(i)=sci;

end

% plot
x=linspace(min(rep1minusmus),max(rep1minusmus));
y=x;

allmus1=[rep1minusmus' rep1plusmus'];
allsigs1=[rep1minussigs' rep1plussigs'];
allmus2=[rep2minusmus' rep2plusmus'];
allsigs2=[rep2minussigs' rep2plussigs'];
mdl = fitlm(allmus1(:,1),allmus2(:,1));
rsq=mdl.Rsquared.Adjusted;

%%
setfig('compare replicates');clf
x=linspace(min(edges),max(edges));
y=x;

hold on
% try
% dscatter(allmus1(:,1),allmus2(:,1),'BINS',[50 50])
% catch
% plot(allmus2(:,1),allmus2(:,2),'o','color',[70 0 0],'facealpha',0.5,'markersize',10)
% end
dscatter(allmus1(:,1), allmus2(:,1))
mdl = fitlm(allmus1(:,1),allmus2(:,1));
rsq=mdl.Rsquared.Adjusted;

% plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
plot(x,x,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('2-')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,numseqs);
text(2,5.5,t)


%%
setfig('switches?');clf
allmus=[allmus1 allmus2];
labelnames={'1-','1+','2-','2+'};

for i=1:4
    for j=1:4
    subplot(4,4,4*(i-1)+j)
    dscatter(allmus(:,i),allmus(:,j))
    xlabel(labelnames{i})
    ylabel(labelnames{j})
    axis([1 6 1 6])
    hold on
    plot(x,x,'k:','linewidth',1.5)
    hold off
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',13)
    end
end
    
    


%% find histogram of distribution of coverage
setfig('coverage dist');clf
[n,c]=hist(tot1minus+tot1plus+tot2minus+tot2plus,length(tot1minus+tot1plus)/250);
area(c,n)
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('coverage')
ylabel('no. of sequences')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)


%% extract sequences
fid=fopen('MFSIV_parsed_uniqcountsseqs.txt');
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


%%
data.origbincounts=mfsIV;
data.scaledbincounts=scaledcounts;
data.seqs=seqs{1}(goodstats);
data.mus=allmus;
data.sigma=[allsigs1 allsigs2];
data.rep1ind=rep1minus;
data.rep2ind=rep2minus;
data.rep3ind=rep1plus;
data.rep4ind=rep2plus;
data.combinedcounts=data.scaledbincounts(:,data.rep1ind)+data.scaledbincounts(:,data.rep2ind);
data.combinedmus=[];
data.combinedsigma=[];
for i=1:length(data.combinedcounts(:,1))
    edges=1:6;
    [m,s]=normfit(edges,[],[],(data.combinedcounts(i,:)));
    data.combinedmus(i)=m;
    data.combinedsigma(i)=s;
end
    
theoapt='ATACCAGCATCGTCTTGATGCCCTTGGAAG';
FAapt='TGCTTGGTACGTTATATTCAG';

data.VYBmus=0.6566*allmus-2.55915;
%% pull out sequences that contain theophylline aptamer

s=regexp(data.seqs,theoapt);
hastheo=~cellfun('isempty',s);
theomus=data.mus(hastheo,:);
theoseqs=data.seqs(hastheo);
theo.mus=theomus;
theo.sigma=data.sigma(hastheo,:);
theo.seqs=theoseqs;
theo.combinedmus=data.combinedmus(hastheo);
theo.combinedsigma=data.combinedsigma(hastheo);
theo.VYBmus=data.VYBmus(hastheo,:);
theo.sigma=data.sigma(hastheo,:);
theo.origbincounts=data.origbincounts(hastheo,:);
theo.scaledbincounts=data.scaledbincounts(hastheo,:);
theo.rep1ind=data.rep1ind;
theo.rep2ind=data.rep2ind;
theo.rep3ind=data.rep3ind;
theo.rep4ind=data.rep4ind;

save('theo.mat','theo');
setfig('theo');clf
labelnames={'1-','Lig+','2-','Lig++'};
% 
% for i=1:4
%     for j=1:4
%     subplot(4,4,4*(i-1)+j)
%     dscatter(theomus(:,i),theomus(:,j))
%     xlabel(labelnames{i})
%     ylabel(labelnames{j})
%     axis([1 6 1 6])
%     hold on
%     plot(x,x,'k:','linewidth',1.5)
%     set(gca,'linewidth',1.5)
%     set(gca,'fontsize',13)
%     title('theo')
%     end
% end

% spiked in controls
sTRSV='GCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGC';
g814='GCTGTCACCGGATGTTTCCGGTCTGATGAGTCCATAAGGACGAAACAGC';
g833='GCTGTCACCGGATACTTCCGGTCTGATGAGTCCCAGAGGACGAAACAGC';
g862='GCTGTCACCGGATGCATCCGGTCTGATGAGTCCCGCGGGACGAAACAGC';
L2b9a1='GCTGTCACCGGAATCAAGGTCCGGTCTGATGAGTCCGTTGTCCAATACCAGCATCGTCTTGATGCCCTTGGCAGTGGATGGGGACGGAGGACGAAACAGC';
Rs2='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAAAGGACGAAACAGC';
Theo1041='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAGAGGACGAAACAGC';
sTRSVctl='GCTGTCACCGGATGTGCTTTCCGGTACGTGAGGTCCGTGAGGACGAAACAGC';

theolibcontrols={sTRSV,g814,g833,g862,L2b9a1,Theo1041,Rs2};
ctrlnames={'sTRSV','g814','g833','g862','L2b9a1','Theo1041','Rs2'};
controlsMFSIV=zeros(length(theolibcontrols),4);

x=linspace(min(min(data.VYBmus)),max(max(data.VYBmus)));
setfig('theo');clf
subplot(1,3,1)
hold on
dscatter(theo.VYBmus(:,1),theo.VYBmus(:,3))
plot(x,x,'k:','linewidth',1.5)

subplot(1,3,2)
hold on
dscatter(theo.VYBmus(:,1),theo.VYBmus(:,2))
plot(x,x,'k:','linewidth',1.5)

subplot(1,3,3)
hold on
dscatter(theo.VYBmus(:,1),theo.VYBmus(:,4))
plot(x,x,'k:','linewidth',1.5)

    
for k=1:length(theolibcontrols)
    s=regexp(data.seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    ctrlmus=data.VYBmus(ctrl,:);
    
    
    setfig('theo')
    subplot(1,3,1)
    plot(ctrlmus(:,1),ctrlmus(:,3),'o','linewidth',2)
    xlabel('Replicate 1 no theophylline')
    ylabel('Replicate 2 no theophylline')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    axis([-1 1.3 -1 1.3])
    
    subplot(1,3,2)
    hold on
    plot(ctrlmus(:,1),ctrlmus(:,2),'o','linewidth',2)
    xlabel('No theophylline')
    ylabel('1 mM theophylline')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    axis([-1 1.3 -1 1.3])
    subplot(1,3,3)
    hold on
    plot(ctrlmus(:,1),ctrlmus(:,4),'o','linewidth',2)
    xlabel('No theophylline')
    ylabel('5 mM theophylline') 
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)
    axis([-1 1.3 -1 1.3])
    controlsMFSIV(k,:)=ctrlmus(1,:);
    

end

legend('Library','1:1',ctrlnames{:},'Location','Best')

%% pull out sequences that contain FA aptamer

s=regexp(data.seqs,FAapt);
hasFA=~cellfun('isempty',s);
FAmus=data.mus(hasFA,:);
FAseqs=data.seqs(hasFA);
FA.mus=FAmus;
FA.sigma=data.sigma(hasFA,:);

FA.seqs=FAseqs;
FA.VYBmus=data.VYBmus(hasFA,:);
FA.sigma=data.sigma(hasFA,:);
FA.combinedmus=data.combinedmus(hasFA);
FA.combinedsigma=data.combinedsigma(hasFA);
FA.origbincounts=data.origbincounts(hasFA,:);
FA.scaledbincounts=data.scaledbincounts(hasFA,:);

FA.rep1ind=data.rep1ind;
FA.rep2ind=data.rep2ind;
FA.rep3ind=data.rep3ind;
FA.rep4ind=data.rep4ind;
save('FA.mat','FA');


setfig('FA');clf

% for i=1:4
%     for j=1:4
%     subplot(4,4,4*(i-1)+j)
%     dscatter(FAmus(:,i),FAmus(:,j))
%     xlabel(labelnames{i})
%     ylabel(labelnames{j})
%     axis([1 6 1 6])
%     hold on
%     plot(x,x,'k:','linewidth',1.5)
%     hold off
%     set(gca,'linewidth',1.5)
%     set(gca,'fontsize',13)
%     title('FA')
%     end
% end



% spiked in controls
FAtert14='GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCTTGGAGAGACGAAACAGC';

s=regexp(data.seqs,FAtert14);
ctrl=~cellfun('isempty',s);
ctrlmus=data.mus(ctrl,:);

setfig('FA')
subplot(1,3,1)
hold on
dscatter(FAmus(:,1),FAmus(:,3))
plot(ctrlmus(:,1),ctrlmus(:,3),'o','linewidth',2)
plot(linspace(min(FAmus(:,1)),max(FAmus(:,1))),linspace(min(FAmus(:,1)),max(FAmus(:,1))),'k:','linewidth',1.5)
xlabel('Replicate 1 no folinic acid')
ylabel('Replicate 2 no folinic acid')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
    axis([2 5.5 2 5.5])
subplot(1,3,2)
hold on
dscatter(FAmus(:,1),FAmus(:,2))
plot(ctrlmus(:,1),ctrlmus(:,2),'o','linewidth',2)
plot(linspace(min(FAmus(:,1)),max(FAmus(:,1))),linspace(min(FAmus(:,1)),max(FAmus(:,1))),'k:','linewidth',1.5)
xlabel('No folinic acid')
ylabel('2 mM R,S folinic acid')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
    axis([2 5.5 2 5.5])
subplot(1,3,3)
hold on
dscatter(FAmus(:,1),FAmus(:,4))
plot(ctrlmus(:,1),ctrlmus(:,4),'o','linewidth',2)
plot(linspace(min(FAmus(:,1)),max(FAmus(:,1))),linspace(min(FAmus(:,1)),max(FAmus(:,1))),'k:','linewidth',1.5)
xlabel('No folinic acid')
ylabel('6 mM R,S folinic acid') 
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
    axis([2 5.5 2 5.5])

%% pull out sequences that contains neither aptamer

rbzmus=data.mus(find(~hasFA .* ~hastheo),:);
rbz.mus=rbzmus;
rbz.combinedmus=data.combinedmus(find(~hasFA .* ~hastheo));
rbz.combinedsigma=data.combinedsigma(find(~hasFA .* ~hastheo));
rbz.seqs=data.seqs(find(~hasFA .* ~hastheo),:);
rbz.VYBmus=data.VYBmus(find(~hasFA .* ~hastheo),:);
rbz.sigma=data.sigma(find(~hasFA .* ~hastheo),:);
rbz.origbincounts=data.origbincounts(find(~hasFA .* ~hastheo),:);
rbz.scaledbincounts=data.scaledbincounts(find(~hasFA .* ~hastheo),:);

rbz.rep1ind=data.rep1ind;
rbz.rep2ind=data.rep2ind;
rbz.rep3ind=data.rep3ind;
rbz.rep4ind=data.rep4ind;



save('rbz.mat','rbz');

setfig('Rbz');clf
subplot(1,3,1)
hold on
dscatter(rbzmus(:,1),rbzmus(:,3))
plot(x,x,'k:','linewidth',1.5)
xlabel(labelnames{1})
ylabel(labelnames{3})
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
    axis([2 5.5 2 5.5])

subplot(1,3,2)
hold on
dscatter(rbzmus(:,1),rbzmus(:,2))
plot(x,x,'k:','linewidth',1.5)
xlabel(labelnames{1})
ylabel(labelnames{2})
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
    axis([2 5.5 2 5.5])

subplot(1,3,3)
hold on
dscatter(rbzmus(:,1),rbzmus(:,4))
plot(x,x,'k:','linewidth',1.5)
xlabel(labelnames{1})
ylabel(labelnames{4}) 
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
    axis([2 5.5 2 5.5])

%%
setfig('make figure MFSIV');clf
hold on
dscatter(data.VYBmus(:,1),data.VYBmus(:,3),'BINS',[50 50])
% plot(theo.VYBmus(:,1),theo.VYBmus(:,2),'.','markersize',10)
% plot(FA.VYBmus(:,1),FA.VYBmus(:,2),'.','markersize',10)

% ctrlmus=zeros(length(theolibcontrols),2);
% for k=1:length(theolibcontrols)
%     s=regexp(alldata.seqs,theolibcontrols{k});
%     ctrl=~cellfun('isempty',s);
%     ctrlmus(k,:)=alldata.VYBmus(ctrl,:);
%     plot(ctrlmus(k,1),ctrlmus(k,2),'o','linewidth',2,'markersize',10)
% 
% end
% legend({'Ribozyme library','Theophylline library','Folinic acid library',ctrlnames{:}},'location','best')
% x=linspace(min(alldata.VYBmus(:,1)),max(alldata.VYBmus(:,2)));
% % plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
plot(x,x,'k:','linewidth',1.5)
mdl = fitlm(rbz.VYBmus(:,1),rbz.VYBmus(:,3));
rsq=mdl.Rsquared.Adjusted;

% title('Mammalian')
xlabel('Replicate 1 log10(\mu)')
ylabel('Replicate 2 log10(\mu)')
set(gca,'linewidth',1.5)
set(gca,'fontsize',18)
set(gca,'XTick',-2:0.5:1.5)
t=sprintf('Corr. = %0.2f',sqrt(rsq));
text(-1.5,0.8,t,'fontsize',18)



%% END 






