addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq/
clear

d=load('~/Documents/CS273/rbz.mat');
allrbz=d.rbz;
x=linspace(min(min(allrbz.VYBmus)),max(max(allrbz.VYBmus)));
labelnames={'1','2'};
numseqs=10;
allrbz(1).motifname='all ribozymes';

data=findSEMbyBootstrp(allrbz,10,1:12);

%%
myttest=struct;
myttest.H=[];
myttest.P=[];
[H,P]=ttest2(data.bootcountsR1',data.bootcountsR2',1e-2,'both','unequal');

%% compare sigmas and sem
setfig('sem');clf
subplot(2,1,1)
[n,c]=hist(allrbz.semC,100);
area(c,n)
title('SEM')

subplot(2,1,2)
[n,c]=hist(allrbz.combinedsigma,100);
hold on
area(c,n)
title('STD')







%% END








