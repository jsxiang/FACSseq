addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq/
clear 
c=load('common.mat');
common=c.common;
numseqs=125;

for k=1:length(common)
%     common(k).origbincounts(:,common(k).rep1ind)=common(k).origbincounts(:,common(k).rep1ind)+common(k).origbincounts(:,common(k).rep1ind+2);
%     common(k).origbincounts(:,common(k).rep1ind)=common(k).origbincounts(:,common(k).rep1ind+2);
    common(k).num1=sum(common(k).origbincounts(:,common(k).rep1ind),2);
    common(k).num2=sum(common(k).origbincounts(:,common(k).rep2ind),2);
    common(k).goodnum=common(k).num1>(3*numseqs) & common(k).num2>numseqs;
    commondata(k)=findSEMbyBootstrp(common(k),50);
    
    
    common1n2(k)=common(k);
    common1n2(k).origbincounts(:,common(k).rep1ind)=common(k).origbincounts(:,common(k).rep1ind)+common(k).origbincounts(:,common(k).rep1ind+1);
    commondata1n2(k)=findSEMbyBootstrp(common1n2(k),50);
    commondata1n2(k).num1=sum(commondata1n2(k).origbincounts(:,commondata1n2(k).rep1ind),2);
    commondata1n2(k).goodnum=commondata1n2(k).num1>(3*numseqs) & commondata1n2(k).num2>numseqs;
    
    common2(k)=common(k);
    common2(k).origbincounts(:,common(k).rep1ind)=common(k).origbincounts(:,common(k).rep1ind+1);
    commondata2(k)=findSEMbyBootstrp(common2(k),50);
    commondata2(k).num1=sum(commondata2(k).origbincounts(:,commondata2(k).rep1ind),2);
    commondata2(k).goodnum=commondata2(k).num1>(3*numseqs) & commondata2(k).num2>numseqs;
    
    commonplus(k)=common(k);
    if k==2
        commonplus(k).rep1ind=common(k).rep1ind+3;
    else 
        commonplus(k).rep1ind=common(k).rep1ind+2;
    end
    commonplus(k).rep2ind=common(k).rep2ind+2;
    commondataplus(k)=findSEMbyBootstrp(commonplus(k),50);
    
end
commondata1=commondata;
%%
setfig('replicates?');clf


subplot(1,3,1)
for k=1:length(commondata1)
    hold on
    plot(mean(commondata1(k).bootcountsR1,2),mean(commondata2(k).bootcountsR1,2),'.')
end

subplot(1,3,2)
for k=1:length(commondata1)
    hold on
    plot(mean(commondata1(k).bootcountsR1,2),mean(commondata1n2(k).bootcountsR1,2),'.')
end

subplot(1,3,3)
for k=1:length(commondata2)
    hold on
    plot(mean(commondata2(k).bootcountsR1,2),mean(commondata1n2(k).bootcountsR1,2),'.')
end




%% 
setfig('counts');clf
hold on
for k=1:length(common)
    plot(common(k).num1, common(k).num2,'.','markersize',12)   
    median(common(k).num1)
    median(common(k).num2)
end    
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'Fontsize',18)
set(gca,'Linewidth',1.5)
xlabel('Mammalian FACSseq count')
ylabel('Yeast FACSseq count')


%%
setfig('deviants');clf
setfig('switches');clf
lib={'Theophylline','Folinic acid','Ribozymes'};

for k=1:length(commondata1n2)
    commondata1n2(k).num1=sum(commondata1n2(k).origbincounts(:,commondata1n2(k).rep1ind),2);
    commondata1n2(k).goodnum=commondata1n2(k).num1>(3*numseqs) & commondata1n2(k).num2>numseqs;
end
commondata=commondata1n2;


theocoeff1=2.20;
theocoeff2=-0.12;


sTRSV='GCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGC';
g814='GCTGTCACCGGATGTTTCCGGTCTGATGAGTCCATAAGGACGAAACAGC';
g833='GCTGTCACCGGATACTTCCGGTCTGATGAGTCCCAGAGGACGAAACAGC';
g862='GCTGTCACCGGATGCATCCGGTCTGATGAGTCCCGCGGGACGAAACAGC';
L2b9a1='GCTGTCACCGGAATCAAGGTCCGGTCTGATGAGTCCGTTGTCCAATACCAGCATCGTCTTGATGCCCTTGGCAGTGGATGGGGACGGAGGACGAAACAGC';
Rs2='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAAAGGACGAAACAGC';
Theo1041='GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAGAGGACGAAACAGC';
sTRSVctl='GCTGTCACCGGATGTGCTTTCCGGTACGTGAGGTCCGTGAGGACGAAACAGC';

theolibcontrols={sTRSV,g814,g833,L2b9a1,Theo1041,Rs2};
ctrlnames={'sTRSV','g814','g833','L2b9a1','Theo1041','Rs2'};
controlsmus=zeros(length(theolibcontrols),4);
controlssigmas=zeros(length(theolibcontrols),4);

for k=1:length(commondata)
    setfig('deviants')
    commondata(k).rescaledR1=0.6566*commondata(k).bootcountsR1-2.55915;
    commondata(k).rescaledR2=(commondata(k).bootcountsR2-6.7288)/1.9938;
    commondataplus(k).rescaledR1=0.6566*commondataplus(k).bootcountsR1-2.55915;
    commondataplus(k).rescaledR2=(commondataplus(k).bootcountsR2-6.7288)/1.9938;
    
    [H,P]=ttest2(commondata(k).rescaledR1',(commondata(k).rescaledR2'-theocoeff2)/theocoeff1,1e-60,'both','unequal');
    commondata(k).ttestH=H;
    
    [H,P]=ttest2(commondataplus(k).rescaledR1',(commondataplus(k).rescaledR2'-theocoeff2)/theocoeff1,1e-60,'both','unequal');
    commondataplus(k).ttestH=H;
    
    subplot(1,3,k)
    hold on
    xmu=mean(commondata(k).rescaledR1(commondata(k).goodnum,:),2);
    ymu=mean(commondata(k).rescaledR2(commondata(k).goodnum,:),2);
    xerr=std(commondata(k).rescaledR1(commondata(k).goodnum,:),[],2);
    yerr=std(commondata(k).rescaledR2(commondata(k).goodnum,:),[],2);
    
    xmuplus=mean(commondataplus(k).rescaledR1(commondata(k).goodnum,:),2);
    ymuplus=mean(commondataplus(k).rescaledR2(commondata(k).goodnum,:),2);
    xerrplus=std(commondataplus(k).rescaledR1(commondata(k).goodnum,:),[],2);
    yerrplus=std(commondataplus(k).rescaledR2(commondata(k).goodnum,:),[],2);
    
    xdevplus=mean(commondataplus(k).rescaledR1(find(commondataplus(k).ttestH .* commondataplus(k).goodnum'),:),2);
    ydevplus=mean(commondataplus(k).rescaledR2(find(commondataplus(k).ttestH .* commondataplus(k).goodnum'),:),2);
    
    xdev=mean(commondata(k).rescaledR1(find(commondata(k).ttestH .* commondata(k).goodnum'),:),2);
    ydev=mean(commondata(k).rescaledR2(find(commondata(k).ttestH .* commondata(k).goodnum'),:),2);
    
    h=herrorbar(xmu,ymu,xerr,xerr,'.',[0.6 0.6 0.6])
    eh=errorbar(xmu,ymu,yerr,yerr,'.','color',[0.6 0.6 0.6],'MarkerSize',12)

    herrorbar(xmuplus,ymuplus,xerrplus,xerrplus,'.',[0.9 0.6 0.6])
    errorbar(xmuplus,ymuplus,yerrplus,yerrplus,'.','color',[0.9 0.6 0.6])

    plot(xdev,ydev,'.','markersize',15,'color',[0.1 0.1 0.1])
    plot(xdevplus,ydevplus,'.','markersize',15,'color',[0.6 0.2 0.2])
    
    x=linspace(-.45,0.9);
    plot(x,x*theocoeff1+theocoeff2,':','linewidth',2,'color',[0.6 0.6 0.6])
   
    if k==3
    mdl = fitlm([xmu' xmuplus'],[ymu' ymuplus'],'weights',1./([xerr' xerrplus'].*[yerr' yerrplus']));
    else
    mdl = fitlm([xmu'],[ymu'],'weights',1./(xerr.*yerr));
    end
    
    rsq=mdl.Rsquared.Adjusted;
    plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',2)
    t=sprintf('R^2 = %0.2f\ny = %0.2fx + %0.2f\nn = %0.0f',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1),sum(commondata(k).goodnum));
    if k==3
    text(0.2,-0.55,t,'fontsize',18)
    else
    text(0.3,-0.55,t,'fontsize',18)
    end
    title(lib{k})
    
    xlabel('mammalian log10(mCherry/BFP)')
    ylabel('yeast log10(GFP/mCherry)') 
    set(gca,'linewidth',2)
    set(gca,'fontsize',18)
    axis([-1 1 -1 2.5])
    
    for i=1:length(theolibcontrols)

    s=regexp(commondata(k).seqs,theolibcontrols{i});
    ctrl=~cellfun('isempty',s);
    ctrlmus=[mean(commondata(k).rescaledR1(ctrl,:),2) mean(commondataplus(k).rescaledR1(ctrl,:),2) mean(commondata(k).rescaledR2(ctrl,:),2) mean(commondataplus(k).rescaledR2(ctrl,:),2) ];
    ctrlsigs=[std(commondata(k).rescaledR1(ctrl,:),[],2) std(commondataplus(k).rescaledR1(ctrl,:),[],2) std(commondata(k).rescaledR2(ctrl,:),[],2) std(commondataplus(k).rescaledR2(ctrl,:),[],2) ];
    try
    controlsmus(i,:)=ctrlmus;
    controlssigmas(i,:)=ctrlsigs;
    end
    try
    plot(ctrlmus(:,1),ctrlmus(:,3),'o','linewidth',2,'markersize',10,'color',[0.8 0.4 0.1])
%     plot(ctrlmus(:,2),ctrlmus(:,4),'o','linewidth',2,'markersize',20,'color',[0.8 0.4 0.1])
    end
    end
    
    %%%%%%% switches %%%%%%%
    setfig('switches')
    foldthresh=1;
    
    
    one2onecoeff1=-0.0556;
    one2onecoeff2=1.0595;
    
    [S,P]=ttest2(one2onecoeff1+one2onecoeff2*commondata(k).rescaledR1',commondataplus(k).rescaledR1',1e-50/length(commondata(k).seqs),'both','unequal');
    commondata(k).switches=S .* commondata(k).goodnum';
    
        
    foldchange=10.^(mean(commondataplus(k).rescaledR1,2)-mean(commondata(k).rescaledR1,2));
    commondata(k).goodswitch=find((foldchange>foldthresh)' .* commondata(k).switches);
    
    xswitch=mean(commondata(k).rescaledR1((commondata(k).goodswitch),:),2);
    yswitch=mean(commondataplus(k).rescaledR1((commondata(k).goodswitch),:),2);
    xswitcherr=std(commondata(k).rescaledR1((commondata(k).goodswitch),:),[],2);
    yswitcherr=std(commondataplus(k).rescaledR1((commondata(k).goodswitch),:),[],2);

    
    subplot(1,3,k)
    hold on
    herrorbar(xmu,xmuplus,xerr,xerr,'.',[0.6 0.6 0.6])
    errorbar(xmu,xmuplus,xerrplus,xerrplus,'.','color',[0.6 0.6 0.6],'markersize',15)

    herrorbar(xswitch,yswitch,xswitcherr,xswitcherr,'.',[0.9 0.3 0.3])
    errorbar(xswitch,yswitch,yswitcherr,yswitcherr,'.','color',[0.9 0.3 0.3],'markersize',15)
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',18)

    plot(x,one2onecoeff1+x*one2onecoeff2,'k:','linewidth',2)
    
    if k==3
    mdl = fitlm([xmu'],[xmuplus'],'weights',1./(xerr.*xerr));
    rsq=mdl.Rsquared.Adjusted;
    plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',2)
    t=sprintf('R^2 = %0.2f\ny = %0.2fx + %0.2f\n',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
    text(0.15,-0.35,t,'fontsize',18)
    end
    axis([-0.8 1 -0.5 1])
    title(lib{k})
    xlabel('log10(mCherry/BFP) no ligand')
    ylabel('log10(mCherry/BFP) with ligand')
    
    if k==1
    f=fopen('switches_yeastmamm.fasta','w');
    else
    f=fopen('switches_yeastmamm.fasta','a');
    end
    for i=1:length(commondata(k).goodswitch)
        fprintf(f,'>%s\n%s\n',lib{k},commondata(k).seqs{commondata(k).goodswitch(i)});
    end
    fclose(f);
    
    if k==1
    f=fopen('gooddata.txt','w');

    else
    f=fopen('gooddata.txt','a');
    end
    fprintf(f,'%s hits\n',lib{k})
    fprintf(f,'m.mu-\tm.SEM-\tm.mu+\tm.SEM+\tm.fold\ty.mu-\ty.SEM-\ty.mu+\ty.SEM+\ty.fold\tsequence\n')
    
    switchdata(k).mamm.muminus=mean(commondata(k).rescaledR1(commondata(k).goodswitch,:),2);
    switchdata(k).mamm.semminus=std(commondata(k).rescaledR1(commondata(k).goodswitch,:),[],2);
    switchdata(k).mamm.muplus=mean(commondataplus(k).rescaledR1(commondata(k).goodswitch,:),2);
    switchdata(k).mamm.semplus=std(commondataplus(k).rescaledR1(commondata(k).goodswitch,:),[],2);
    switchdata(k).mamm.fold=10.^(mean(commondataplus(k).rescaledR1(commondata(k).goodswitch,:),2)-mean(commondata(k).rescaledR1(commondata(k).goodswitch,:),2));
    switchdata(k).mamm.seqs=commondata(k).seqs(commondata(k).goodswitch);
    
    switchdata(k).yeast.muminus=mean(commondata(k).rescaledR2(commondata(k).goodswitch,:),2);
    switchdata(k).yeast.semminus=std(commondata(k).rescaledR2(commondata(k).goodswitch,:),[],2);
    switchdata(k).yeast.muplus=mean(commondataplus(k).rescaledR2(commondata(k).goodswitch,:),2);
    switchdata(k).yeast.semplus=std(commondataplus(k).rescaledR2(commondata(k).goodswitch,:),[],2);
    switchdata(k).yeast.fold=10.^(mean(commondataplus(k).rescaledR2(commondata(k).goodswitch,:),2)-mean(commondata(k).rescaledR2(commondata(k).goodswitch,:),2));
    switchdata(k).yeast.seqs=commondata(k).seqs(commondata(k).goodswitch(i));
    
    gooddata(k).mamm.muminus=mean(commondata(k).rescaledR1(commondata(k).goodnum,:),2);
    gooddata(k).mamm.semminus=std(commondata(k).rescaledR1(commondata(k).goodnum,:),[],2);
    gooddata(k).mamm.muplus=mean(commondataplus(k).rescaledR1(commondata(k).goodnum,:),2);
    gooddata(k).mamm.semplus=std(commondataplus(k).rescaledR1(commondata(k).goodnum,:),[],2);
    gooddata(k).mamm.fold=10.^(mean(commondataplus(k).rescaledR1(commondata(k).goodnum,:),2)-mean(commondata(k).rescaledR1(commondata(k).goodnum,:),2));
    gooddata(k).mamm.seqs=commondata(k).seqs(commondata(k).goodnum);
    
    gooddata(k).yeast.muminus=mean(commondata(k).rescaledR2(commondata(k).goodnum,:),2);
    gooddata(k).yeast.semminus=std(commondata(k).rescaledR2(commondata(k).goodnum,:),[],2);
    gooddata(k).yeast.muplus=mean(commondataplus(k).rescaledR2(commondata(k).goodnum,:),2);
    gooddata(k).yeast.semplus=std(commondataplus(k).rescaledR2(commondata(k).goodnum,:),[],2);
    gooddata(k).yeast.fold=10.^(mean(commondataplus(k).rescaledR2(commondata(k).goodnum,:),2)-mean(commondata(k).rescaledR2(commondata(k).goodnum,:),2));
    gooddata(k).yeast.seqs=commondata(k).seqs(commondata(k).goodswitch(i));
    
    
    for i=1:length(commondata(k).goodswitch)
        fprintf(f,'%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%s\n',...
            gooddata(k).mamm.muminus(i),...
            gooddata(k).mamm.semminus(i),...
            gooddata(k).mamm.muplus(i),...
            gooddata(k).mamm.semplus(i),...
            gooddata(k).mamm.fold(i),...
            gooddata(k).yeast.muminus(i),...
            gooddata(k).yeast.semminus(i),...
            gooddata(k).yeast.muplus(i),...
            gooddata(k).yeast.semplus(i),...
            gooddata(k).yeast.fold(i),...
            gooddata(k).mamm.seqs{i});
    end
    fclose(f);

end

controlsmus(controlsmus==0)=nan;
controlssigmas(controlsmus==0)=nan;



%%
setfig('controls only');clf
hold on
for i=1:length(controlsmus(:,1))
    plot(controlsmus(i,1),controlsmus(i,3),'o','markersize',10)
    
    
end
legend(ctrlnames,'location','northwest')
herrorbar(controlsmus(:,1),controlsmus(:,3),controlssigmas(:,1),controlssigmas(:,1),'.',[0.5 0.5 1])
errorbar(controlsmus(:,1),controlsmus(:,3),controlssigmas(:,3),controlssigmas(:,3),'.','color',[0.5 0.5 1])


mdl = fitlm(controlsmus(:,1),controlsmus(:,3));
rsq=mdl.Rsquared.Adjusted;
x=linspace(-0.5,0.5);
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',2)
t=sprintf('R^2 = %0.2f\ny = %0.2fx + %0.2f',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
text(0.1,-0.4,t,'fontsize',18)
xlabel('mammalian log10(mCherry/BFP)')
ylabel('yeast log10(GFP/mCherry)') 
set(gca,'linewidth',2)
set(gca,'fontsize',18)
xlim([-1 1])


setfig('check coverage mamm');clf
hold on
for k=1:length(commondata)
    subplot(3,1,k)
histogram(commondata(k).num1,1:20:1000)
title(lib{k})
xlabel('read count')
ylabel('no. of sequences')
end

setfig('check coverage yeast');clf
hold on
for k=1:length(commondata)
    subplot(3,1,k)
histogram(commondata(k).num2,1:10:350)
title(lib{k})
xlabel('read count')
ylabel('no. of sequences')
end

%% switches

setfig('compare fold change');clf
for k=1:length(switchdata)
    plot(log10(switchdata(k).mamm.fold),log10(switchdata(k).yeast.fold),'.')
end

%% pick hits
setfig('pick hits');clf
numswitches=4;
numrbz=2;
hitseqs={};

handpickedhits(1).mamm.muminus=nan;
handpickedhits(2).mamm.muminus=nan;
handpickedhits(3).mamm.muminus=nan;
for k=1:length(gooddata)
murange=max(gooddata(k).mamm.muminus)-min(gooddata(k).mamm.muminus);
if k==3
mucand=linspace(-0.2062,0.98*murange+min(gooddata(k).mamm.muminus),8);
else
mucand=linspace(0.5*murange+min(gooddata(k).mamm.muminus),0.98*murange+min(gooddata(k).mamm.muminus),numrbz);
end

mucand=[mucand handpickedhits(k).mamm.muminus];

mucandind=[];
for i=1:length(mucand)
    if ~isnan(mucand(i))
    candind=find(abs(mucand(i)-gooddata(k).mamm.muminus)==min(abs(mucand(i)-gooddata(k).mamm.muminus)));    
    mucandind(end+1)=candind(1);
    end
end


hold on
plot(gooddata(k).mamm.muminus(mucandind),gooddata(k).mamm.muplus(mucandind),'o','MarkerSize',10,'linewidth',2)
hitseqs={hitseqs{:} gooddata(k).mamm.seqs{mucandind}};


try
[maxfolds,mfI]=sort(switchdata(k).mamm.fold);
maxfolds=maxfolds(end:-1:end-numswitches);
mfI=mfI(end:-1:end-numswitches);

hitseqs={hitseqs{:} switchdata(k).mamm.seqs{mfI}};
end
end
% setfig('where are the mus');clf
% bar([rep1mus(goodagreement(mucandind));rep2mus(goodagreement(mucandind))]')
% seqcand=gooddata.seqs(goodagreement(mucandind))

uniquehits=unique(hitseqs);



f=fopen('pickedhits.txt','w');
orderedhits={'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAGAAAGGACGAAACAGC',
'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAGGAAGGACGAAACAGC',
'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCATAAGGACGAAACAGC',
'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCTGAAGGACGAAACAGC',
'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCTTGAGGACGAAACAGC',
'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCTGGTAGGACGAAACAGC',
'GCTGTCACCGGGTGCTTGGTACGTTATATTCAGCCCGGTCTGATGAGTCCCGATAGGACGAAACAGC',
'GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCCAAACGGGACGAAACAGC',
'GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCCAAGAGGGACGAAACAGC',
'GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCCGTGAGGGACGAAACAGC',
'GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCCTTAGGGGACGAAACAGC',
'GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCCTTTGGGGACGAAACAGC',
'GCTGTCACCGGTGCTTGGTACGTTATATTCAGCCGGTCTGATGAGTCCCCCCGACGAAACAGC',
'GCTGTCACCGGAGTACCTGTCGACTGTGTGGACAAACATACATCCGGTCTGATGAGTCCTGAAATGGACGAAACAGC',
'GCTGTCACCGGATGTGGAAGTCCGGTCTGATGAGTCCCAGACTACCATACATAAGAGAAACACGCCAGGACGAAACAGC',
'GCTGTCACCGGATACAAAAGTCCGGTCTGATGAGTCCCGGAGACATGCATTCCCATTTATCCTTTTTAGGACGAAACAGC',
'GCTGTCACCGGATGTTCGATATCCCATCGGTGTCGCCTCTAGTCCGGTCTGATGAGTCCCCAAACTAGGACGAAACAGC',
'GCTGTCACCGGATCGCAGCCCTGGGGCAGGCCCATCCGTCAGTCCGGTCTGATGAGTCCGAGGAAGGACGAAACAGC',
'GCTGTCACCGGATTTCGCGTCCGGTCTGATGAGTCCCGCGACGATCCCATTTCAAGTGGTGTAAAAGGACGAAACAGC'};
for k=1:length(gooddata)
fprintf(f,'%s hits\n',lib{k})
fprintf(f,'m.mu-\tm.SEM-\tm.mu+\tm.SEM+\tm.fold\ty.mu-\ty.SEM-\ty.mu+\ty.SEM+\ty.fold\tsequence\n')
   for i=1:length(orderedhits)

    s=regexp(gooddata(k).mamm.seqs,orderedhits{i});
    phits=find(~cellfun('isempty',s));
  
        fprintf(f,'%1.4f\t%1.4f\t%1.4f\t%1.4f\t%2.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%2.4f\t%s\n',...
            gooddata(k).mamm.muminus(phits),...
            gooddata(k).mamm.semminus(phits),...
            gooddata(k).mamm.muplus(phits),...
            gooddata(k).mamm.semplus(phits),...
            gooddata(k).mamm.fold(phits),...
            gooddata(k).yeast.muminus(phits),...
            gooddata(k).yeast.semminus(phits),...
            gooddata(k).yeast.muplus(phits),...
            gooddata(k).yeast.semplus(phits),...
            gooddata(k).yeast.fold(phits),...
            gooddata(k).mamm.seqs{phits});
    end
end
fclose(f);



setfig('deviants')
for k=1:length(commondata)

    for i=1:length(orderedhits)
    setfig('deviants')
    subplot(1,3,k)
    s=regexp(commondata(k).seqs,orderedhits{i});
    hits=~cellfun('isempty',s);
    plot(mean(commondata(k).rescaledR1(hits,:),2),mean(commondata(k).rescaledR2(hits,:),2),'o','color',[0.7 0.1 0.9],'linewidth',2,'markersize',10)
    

    setfig('switches')
    subplot(1,3,k)
    plot(mean(commondata(k).rescaledR1(hits,:),2),mean(commondataplus(k).rescaledR1(hits,:),2),'o','color',[0.7 0.1 0.9],'linewidth',2,'markersize',10)
    end
end

%% check controls
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

hold on
for i=1:length(common) 
        setfig(lib{i})
    for k=1:length(theolibcontrols)

    s=regexp(common(i).seqs,theolibcontrols{k});
    ctrl=~cellfun('isempty',s);
    
    ctrlmus=[commondata(i).mVYBmus1(ctrl)' common(i).mVYBmus2(ctrl)' common(i).yVYBmus1(ctrl)' common(i).yVYBmus2(ctrl)'];
    try
        controlsMFSIV(k,:)=ctrlmus;
    catch
    end
    hold on
%     if i==3
%         plot(ctrlmus(:,1),mean(ctrlmus(:,[3 4]),2),'o','linewidth',2,'markersize',16,'color',[0.8 0.4 0.1])
%     else
    plot(ctrlmus(:,1),ctrlmus(:,3),'o','linewidth',2,'markersize',20,'color',[0.8 0.4 0.1])
    plot(ctrlmus(:,2),ctrlmus(:,4),'o','linewidth',2,'markersize',20,'color',[0.8 0.4 0.1])
    legend('no ligand','with ligand','location','best')
%     end
    end
end

controlsMFSIV(controlsMFSIV==0)=nan;
%%
setfig('controls yeast vs mammalian');clf

% plot(controlsMFSIV(:,1),[mean(controlsMFSIV(1:3,3:4),2); controlsMFSIV(4:7,3)],'.','MarkerSize',15)
plot(controlsMFSIV(:,1),controlsMFSIV(:,3),'.','MarkerSize',20)
mdl = fitlm(controlsMFSIV(:,1),controlsMFSIV(:,3));
rsq=mdl.Rsquared.Adjusted;
xlabel('mammalian cells')
ylabel('yeast')
title('Spiked in controls')
set(gca,'linewidth',2)
set(gca,'fontsize',28)

x=linspace(min(controlsMFSIV(:,1)),max(controlsMFSIV(:,1)));
hold on
plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',2)
t=sprintf('R^2 = %0.2f\n y = %0.2fx + %0.2f',rsq,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
text(-0.25,0.1,t,'fontsize',28)
% legend('Library','1:1',ctrlnames{:},'Location','Best')
% 
% 

%%
f=fopen('seqs_yeastmamm.fasta','w');
fprintf(f,'>seq\n%s\n',commondata(1).seqs{commondata});
fclose(f);
data=findMotif(common(1),{'(GCTGTC[A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G])([A|C|T|G]*)([A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G]CTGATGAG[A|C|T|G][A|C|T|G][A|C|T|G])([A|C|T|G]*)([A|C|T|G][A|C|T|G][A|C|T|G]CGAAACAGC)'},{'rbz'});
data=findSeqnum(data);
l1l2seqs=[];
l1l2seqsSIG=[];
l1l2seqsNOTSIG=[];
for i=1:length(data.loop2seqnum)
    if length(data.loop2seqnum{i})==5 & s(i)
        l1l2seqsSIG=[l1l2seqsSIG;data.loop2seqnum{(i)}];
    elseif length(data.loop2seqnum{i})==5 & ~s(i)
        l1l2seqsNOTSIG=[l1l2seqsNOTSIG;data.loop2seqnum{(i)}];
    end
end
[entlo,plo]=findEnt(l1l2seqsSIG,{'significantly different'})
[enthi,phi]=findEnt(l1l2seqsNOTSIG,{'not different'})
[dent,dp]=findDeltaEnt(entlo,enthi,plo,phi,{'yeast diff mamm'})
% seqlogo(data.loop2(data.loop2len==5))



%% 
setfig('dummy')
%%

pdfsavefig('all');

