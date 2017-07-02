addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq
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
yfsIVRbz=load('YFSIVc_parsed_uniqcountsonly.txt');

%% plot distribution of barcodes
% count_eachbarcodeIII=sum(yfsIVRbz)/sum(sum(yfsIVRbz));
count_eachbarcodeIII=sum(yfsIVRbz);
cellsorted=[3503518	580763	3886643	549534	2456115	2454267	2699072	2450899	540513	1891628	590651	1879294	50722	301089	55751	292456	2190	21978	2613	20457	565	1284	764	1376];
allBarcodeCounts=[2428432	533453	2634845	563491	1609481	2064703	1926662	1970091	479765	1819480	344653	1568409	35117	281912	47346	279033	5349	23450	9223	20189	2264	7359	4822	18883]; % found using countBarcodes.sh
% cellsorted=cellsorted/sum(cellsorted);
cellsorted=cellsorted;
setfig('barcode distribution');clf
plot(cellsorted,count_eachbarcodeIII,'.','MarkerSize',15)

xlabel('cells sorted')
ylabel('NGS count')


hold on
xval=linspace(min(cellsorted),max(cellsorted));
yval=xval;
plot(xval,yval)
ylow=xval-0.02;
yhigh=xval+0.02;
plot(xval,ylow,xval,yhigh)

barcodestoohigh=find(count_eachbarcodeIII>(cellsorted+0.02));
barcodestoolow=find(count_eachbarcodeIII<(cellsorted-0.02));
scalefactor=cellsorted./count_eachbarcodeIII;
scalefactor=cellsorted./cellsorted;
scaledcounts=yfsIVRbz.*repmat(scalefactor,length(yfsIVRbz(:,1)),1);
save('YFSIV_scaled.mat','scaledcounts')
hold off
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
title('YFSIV')



% %% pick sequences with good statistics
% rep1minus=[1 2 5 6 9 10 13 14 17 18 21 22];
% rep1minus=rep1minus(end:-1:1);
% rep1plus=rep1minus+2;
% tot1minus=sum(yfsIVRbz(:,rep1minus),2);
% tot1plus=sum(yfsIVRbz(:,rep1plus),2);
% 
% rsquared=[];
% num=[];
% for n=10:2:250
%     goodstats1=tot1minus>n & tot1plus>n;
% 
%     goodscaledcounts1=scaledcounts(find(goodstats1),:);
% 
%     %fit data to normdist
% 
%     rep1minusmus=zeros(1,length(goodscaledcounts1(:,1)));
%     rep1plusmus=rep1minusmus;
%     for i=1:length(goodscaledcounts1(:,1))
%         c1=goodscaledcounts1(i,rep1minus);
%         edges=1:12;
%         [m,s]=normfit(edges,[],[],c1);
%         rep1minusmus(i)=m;
% 
%         c2=goodscaledcounts1(i,rep1plus);
%         [m,s]=normfit(edges,[],[],c2);
%         rep1plusmus(i)=m;
%     end
% 
% 
% 
%     allmus=[rep1minusmus' rep1plusmus'];
% 
%     mdl = fitlm(allmus(:,1),allmus(:,2));
%     rsquared(end+1)=mdl.Rsquared.Adjusted;
%     num(end+1)=sum(goodstats1);
% end
% 
% x=linspace(min(edges),max(edges));
% y=x;
% 
% setfig('n versus rsquared');clf
% hold on
% plot(10:2:250,rsquared,'linewidth',1.5)
% plot(10:2:250,num/length(goodstats1),'linewidth',1.5)
% xlabel('n')
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)
% legend('r^2','fraction of seqs','location','best')
% title('YFSIV')

%% Choose a numseq, plot data for YFSIV

numseqs=20;
rep1minus=[2 1 6 5 10 9 14 13 18 17 22 21];
% rep1minus=[1 5 9 13 17 21];
rep1minus=rep1minus(end:-1:1);
rep1plus=rep1minus+2;
tot1minus=sum(yfsIVRbz(:,rep1minus),2);
tot1plus=sum(yfsIVRbz(:,rep1plus),2);

goodstats1=tot1minus>numseqs & tot1plus>numseqs;
goodscaledcounts1=scaledcounts(find(goodstats1),:);
origcounts1=yfsIVRbz(find(goodstats1),:);

%fit data to normdist
rep1minusmus=zeros(1,length(goodscaledcounts1(:,1)));
rep1plusmus=rep1minusmus;
rep1minussigs=rep1minusmus;
rep1plussigs=rep1minusmus;
rep1mci=rep1minusmus;
rep2mci=rep1minusmus;
rep1sci=rep1minusmus;
rep2sci=rep1minusmus;
for i=1:length(goodscaledcounts1(:,1))
    c1=goodscaledcounts1(i,rep1minus);
    edges=1:12;
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
end

% plot
x=linspace(min(rep1minusmus),max(rep1minusmus));
y=x;

allmus1=[rep1minusmus' rep1plusmus'];
allsigs1=[rep1minussigs' rep1plussigs'];
allmci=[rep1mci' rep2mci'];

% mdl = fitlm(allmus1(:,1),allmus1(:,2));
% rsq=mdl.Rsquared.Adjusted;
%%
setfig('compare replicates');clf

hold on
plot(allmus1(:,1),allmus1(:,2),'.')

% try
% dscatter(allmus1(:,1),allmus1(:,2),'BINS',[50 50])
% catch
%     plot(allmus1(:,1),allmus1(:,2),'.')
% end
%%
% plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
plot(x,x,'r:','linewidth',1.5)
hold off
xlabel('Rep 1 GFP/mCherry')
ylabel('Rep 2 GFP/mCherry')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
% t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,numseqs);
% text(5,11.5,t)

%% troubleshoot

figure(1);clf;
subplot(1,2,1)
bar(yfsIVRbz(12000,rep1minus))
subplot(1,2,2)
bar(yfsIVRbz(12000,rep1plus))

%%



%% find histogram of distribution of coverage
% setfig('coverage dist');clf
% [n,x]=hist(tot1minus+tot1plus,length(tot1minus+tot1plus)/250);
% area(x,n)
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% xlabel('coverage')
% ylabel('no. of sequences')
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',14)


%% extract sequences
fid=fopen('YFSIVc_parsed_uniqcountsseqs.txt');
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

alldata.origbincounts=yfsIVRbz;
alldata.rep1ind=rep1minus;
alldata.rep2ind=rep1plus;
alldata.scaledbincounts=scaledcounts;
alldata.seqs=seqs{1}(goodstats1);
alldata.mus=allmus1;



% alldata.mci=allmci;
alldata.sigma=allsigs1;
theoapt='ATACCAGCATCGTCTTGATGCCCTTGGAAG';
FAapt='TGCTTGGTACGTTATATTCAG';


%% pull out sequences that contain theophylline aptamer

s=regexp(alldata.seqs,theoapt);
hastheo=~cellfun('isempty',s);
theomus=alldata.mus(hastheo,:);
theoseqs=alldata.seqs(hastheo);
theo.mus=theomus;
theo.seqs=theoseqs;
theo.sigma=alldata.sigma(hastheo,:);
theo.origbincounts=alldata.origbincounts(hastheo,:);
theo.scaledbincounts=alldata.scaledbincounts(hastheo,:);
theo.rep1ind=alldata.rep1ind;
theo.rep2ind=alldata.rep2ind;
save('theo.mat','theo');

labelnames={'1','2'};


% spiked in controls
sTRSV='GCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGC';
g814='GCTGTCACCGGATGTTTCCGGTCTGATGAGTCCATAAGGACGAAACAGC';
g833='GCTGTCACCGGATACTTCCGGTCTGATGAGTCCCAGAGGACGAAACAGC';
% g862='GCTGTCACCGGATGCATCCGGTCTGATGAGTCCCGCGGGACGAAACAGC';
L2b9a1='GCTGTCACCGGAATCAAGGTCCGGTCTGATGAGTCCGTTGTCCAATACCAGCATCGTCTTGATGCCCTTGGCAGTGGATGGGGACGGAGGACGAAACAGC';
Rs2='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAAAGGACGAAACAGC';
Theo1041='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAGAGGACGAAACAGC';
% sTRSVctl='GCTGTCACCGGATGTGCTTTCCGGTACGTGAGGTCCGTGAGGACGAAACAGC';

theolibcontrols={sTRSV,g814,g833,L2b9a1,Rs2,Theo1041};
ctrlnames={'sTRSV','g814','g833','L2b9a1','Rs2','Theo1041'};

%% pull out sequences that contain FA aptamer

s=regexp(alldata.seqs,FAapt);
hasFA=~cellfun('isempty',s);
FAmus=alldata.mus(hasFA,:);
FAseqs=alldata.seqs(hasFA);
FA.mus=FAmus;
FA.seqs=FAseqs;
FA.sigma=alldata.sigma(hasFA,:);
FA.origbincounts=alldata.origbincounts(hasFA,:);
FA.scaledbincounts=alldata.scaledbincounts(hasFA,:);
FA.rep1ind=alldata.rep1ind;
FA.rep2ind=alldata.rep2ind;
save('FA.mat','FA');


setfig('FA');clf

% spiked in controls
FAtert14='GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCTTGGAGAGACGAAACAGC';
s=regexp(alldata.seqs,FAtert14);
ctrl=~cellfun('isempty',s);
ctrlmus=alldata.mus(ctrl,:);

setfig('FA')
hold on
% dscatter(FAmus(:,1),FAmus(:,2))
plot(FAmus(:,1),FAmus(:,2),'.')
plot(ctrlmus(:,1),ctrlmus(:,2),'o','linewidth',2)
plot(x,x,'k:','linewidth',1.5)
xlabel(labelnames{1})
ylabel(labelnames{2})
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
[t,s]=findSeqsTTest(FA,1,0.05,1);
FA.ttestdiff=[t' s'];

hold on
plot(FA.mus(find(FA.ttestdiff(:,2)),1),FA.mus(find(FA.ttestdiff(:,2)),2),'s')
%%
setfig('theo');clf
hold on
% dscatter(theomus(:,1),theomus(:,2))
plot(theomus(:,1),theomus(:,2),'.')
plot(x,x,'k:','linewidth',1.5)

ctrlmus=zeros(6,2);
for k=1:length(theolibcontrols)
    s=regexp(alldata.seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    ctrlmus(k,:)=alldata.mus(ctrl,:) ;
      
    plot(ctrlmus(k,1),ctrlmus(k,2),'o','linewidth',2)
end
xlabel(labelnames{1})
ylabel(labelnames{2})
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
legend({'Library','1:1',ctrlnames{:}},'location','best')
[t,s]=findSeqsTTest(theo,1,0.05,1);
theo.ttestdiff=[t' s'];

hold on
plot(theo.mus(find(theo.ttestdiff(:,2)),1),theo.mus(find(theo.ttestdiff(:,2)),2),'s')
%% pull out sequences that contains neither aptamer

rbzmus=alldata.mus(find(~hasFA .* ~hastheo),:);
rbz.mus=rbzmus;
rbz.seqs=alldata.seqs(find(~hasFA .* ~hastheo),:);
rbz.sigma=alldata.sigma(find(~hasFA .* ~hastheo),:);
rbz.origbincounts=alldata.origbincounts(find(~hasFA .* ~hastheo),:);
rbz.scaledbincounts=alldata.scaledbincounts(find(~hasFA .* ~hastheo),:);
rbz.rep1ind=alldata.rep1ind;
rbz.rep2ind=alldata.rep2ind;
save('rbz.mat','rbz');

setfig('Rbz');clf

hold on
plot(rbzmus(:,1),rbzmus(:,2),'.','markersize',10,'color',[0.6 0.6 0.6])
mdl = fitlm(rbzmus(:,1),rbzmus(:,2));
rsq=mdl.Rsquared.Adjusted;

xlabel(labelnames{1})
ylabel(labelnames{2})
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,numseqs);
text(4.5,11.5,t)

[t,s]=findSeqsTTest(rbz,1,0.05,1);
rbz.ttestdiff=[t' s'];

hold on
plot(rbz.mus(find(~rbz.ttestdiff(:,2)),1),rbz.mus(find(~rbz.ttestdiff(:,2)),2),'.','MarkerSize',2)

mdl = fitlm(rbz.mus(find(~rbz.ttestdiff(:,2)),1),rbz.mus(find(~rbz.ttestdiff(:,2)),2));
rsq=mdl.Rsquared.Adjusted;
t=sprintf('R^2 = %0.2f\nT-test \\alpha = 0.05/n',rsq);
text(4.5,10.9,t)



%%
controls=[  -8.0000e-01  -6.6300e-01  -7.5300e-02   -5.6500e-01  -2.4300e-01]; %VYB data

x=linspace(min(controls),max(controls));

setfig('calibrate');clf
hold on
plot(controls',ctrlmus([1 2 3 5 6],1),'.','markersize',20)
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
mdl.Coefficients.Estimate
mdl = fitlm(controls',ctrlmus([1 2 3 5 6],1));
rsq=mdl.Rsquared.Adjusted;

xlabel('VYB log10(GFP/mCherry)')
ylabel('NGS bins')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,numseqs);
text(-.7,6.4,t)

alldata.VYBmus=(alldata.mus-mdl.Coefficients.Estimate(1))/mdl.Coefficients.Estimate(2);
theo.VYBmus=(theo.mus-mdl.Coefficients.Estimate(1))/mdl.Coefficients.Estimate(2);
FA.VYBmus=(FA.mus-mdl.Coefficients.Estimate(1))/mdl.Coefficients.Estimate(2);
rbz.VYBmus=(rbz.mus-mdl.Coefficients.Estimate(1))/mdl.Coefficients.Estimate(2);
save('theo.mat','theo');
save('FA.mat','FA');
save('rbz.mat','rbz');
save('alldata.mat','alldata');
%%

setfig('all');clf
theolibcontrols={sTRSV,g814,g833,L2b9a1,Rs2,Theo1041};
ctrlnames={'sTRSV','g814','g833','L2b9a1','Rs2','Theo1041'};
hold on
dscatter(rbz.VYBmus(:,1),rbz.VYBmus(:,2),'BINS',[50 50])
plot(theo.VYBmus(:,1),theo.VYBmus(:,2),'.','markersize',10)
plot(FA.VYBmus(:,1),FA.VYBmus(:,2),'.','markersize',10)

ctrlmus=zeros(length(theolibcontrols),2);
for k=1:length(theolibcontrols)
    s=regexp(alldata.seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    ctrlmus(k,:)=alldata.VYBmus(ctrl,:);
    plot(ctrlmus(k,1),ctrlmus(k,2),'o','linewidth',2,'markersize',10)

end
legend({'Ribozyme library','Theophylline library','Folinic acid library',ctrlnames{:}},'location','best')
x=linspace(min(alldata.VYBmus(:,1)),max(alldata.VYBmus(:,2)));
% plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
plot(x,x,'k:','linewidth',1.5)
mdl = fitlm(rbz.VYBmus(:,1),rbz.VYBmus(:,2));
rsq=mdl.Rsquared.Adjusted;

xlabel('Rep 1')
ylabel('Rep 2')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
t=sprintf('R^2 = %0.2f\nn = %0.0f',rsq,numseqs);
text(-1.3,2,t)




setfig('make figure YFSIVc');clf
hold on
dscatter(rbz.VYBmus(:,1),rbz.VYBmus(:,2),'BINS',[50 50])
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
mdl = fitlm(rbz.VYBmus(:,1),rbz.VYBmus(:,2));
rsq=mdl.Rsquared.Adjusted;

title('Yeast')
xlabel('Rep 1 \mu')
ylabel('Rep 2 \mu')
set(gca,'linewidth',1.5)
set(gca,'fontsize',18)
t=sprintf('R^2 = %0.2f',rsq);
text(-1.3,2,t,'fontsize',18)



%% 
%%

setfig('all good CoV');clf
theolibcontrols={sTRSV,g814,g833,L2b9a1,Rs2,Theo1041};
ctrlnames={'sTRSV','g814','g833','L2b9a1','Rs2','Theo1041'};
hold on
dscatter(rbz.VYBmus(rbz.goodCoV,1),rbz.VYBmus(rbz.goodCoV,2),'BINS',[50 50])
plot(theo.VYBmus(:,1),theo.VYBmus(:,2),'.','markersize',10)
plot(FA.VYBmus(:,1),FA.VYBmus(:,2),'.','markersize',10)

ctrlmus=zeros(length(theolibcontrols),2);
for k=1:length(theolibcontrols)
    s=regexp(alldata.seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    ctrlmus(k,:)=alldata.VYBmus(ctrl,:);
    plot(ctrlmus(k,1),ctrlmus(k,2),'o','linewidth',2,'markersize',10)

end
legend({'Ribozyme library','Theophylline library','Folinic acid library',ctrlnames{:}},'location','best')
x=linspace(min(alldata.VYBmus(:,1)),max(alldata.VYBmus(:,2)));
% plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
plot(x,x,'k:','linewidth',1.5)

%% END 






