function [t1,significantseqs]=findSeqsTTest(data,bonferroni,alpha,findswitch)
%% Look for statistically significant switches with a t test
rep1mus=data.mus(:,1);
rep2mus=data.mus(:,2);
rep1sigma=data.sigma(:,1);
rep2sigma=data.sigma(:,2);
rep1totcount=sum(data.origbincounts(:,data.rep1ind),2);
rep2totcount=sum(data.origbincounts(:,data.rep2ind),2);
n=length(rep1mus);
a=alpha;
if bonferroni
bonferroni=a/n;
end
if findswitch
    allrep1muslin_minus=rep1mus;
else
variancearea = rep1sigma.*rep2sigma;   % suggested by statistics consulting group
% mdlrep1 = fitlm(rep1mus,rep2mus,'weights',1./variancearea);
mdlrep1 = fitlm(rep1mus,rep2mus);
allrep1muslin_minus=mdlrep1.Coefficients.Estimate(1)+mdlrep1.Coefficients.Estimate(2).*rep1mus;
end
t1=[];
for i=1:length(rep2mus)
    t=(allrep1muslin_minus(i)-rep2mus(i))./sqrt(rep1sigma(i)^2/rep1totcount((i))+rep2sigma(i)^2/rep2totcount((i)));
%     t=(rep1mus(i)-rep2mus(i))./sqrt(rep1sigma(i)^2/rep1totcount((i))+rep2sigma(i)^2/rep2totcount((i)));
    t1(end+1)=tpdf(t,rep1totcount((i))+rep2totcount((i))-1);

end
if bonferroni
    significantseqs=t1<bonferroni;
else
    significantseqs=(t1<alpha);
end
end