function data=findSEMbyBootstrp(inputdata,Nboot)
% function takes FACSseq data object and appends new attributes
% semC of combined
% semR1 of replicate 1
% semR2 of replicate 2
% Nboot is the number of bootstrapped samples generated
% edges is the no. of histogram bins
data=inputdata;
data.semC=zeros(length(data.seqs),1);
data.semR1=zeros(length(data.seqs),1);
data.semR2=zeros(length(data.seqs),1);

data.bootcountsC=zeros(length(data.seqs),Nboot);
data.bootcountsR1=zeros(length(data.seqs),Nboot);
data.bootcountsR2=zeros(length(data.seqs),Nboot);

try
    data.scaledbincounts=inputdata.scaledbincounts;
catch
    data.scaledbincounts=inputdata.origbincounts;
end

for i=1:length(data.seqs)
    xdataC=[];
    xdata1=[];
    xdata2=[];
    if length(data.rep1ind)==length(data.rep2ind)
        combinedcounts=data.scaledbincounts(i,data.rep1ind)+data.scaledbincounts(i,data.rep2ind);
        edges=1:length(combinedcounts(1,:));
        for j=edges
            xdataC=[xdataC repmat(edges(j),1,combinedcounts(j))];
        end
        
            
        try
            mboot=bootstrp(Nboot,@normfit,xdataC);
            data.semC(i)=std(mboot);
            data.bootcountsC(i,:)=mboot';
        catch
            fprintf('bootstrap error with sample %0.0f;likely no. bins with counts <2\n',i);
            xdataC
            data.semC(i)=nan;
        end
    end
    rep1counts=data.scaledbincounts(i,data.rep1ind);
    rep2counts=data.scaledbincounts(i,data.rep2ind);

    edges=1:length(data.rep1ind);
    for j=edges
        xdata1=[xdata1 repmat(edges(j),1,rep1counts(j))];
    end
    edges=1:length(data.rep2ind);
    for j=edges
        xdata2=[xdata2 repmat(edges(j),1,rep2counts(j))];
    end


    try
        mboot=bootstrp(Nboot,@normfit,xdata1);
        data.semR1(i)=std(mboot);
        data.bootcountsR1(i,:)=mboot';
    catch
        fprintf('bootstrap error with sample %0.0f;likely no. bins with counts <2\n',i);
        data.semR1(i)=nan;
    end

    try
        mboot=bootstrp(Nboot,@normfit,xdata2);
        data.semR2(i)=std(mboot);
        data.bootcountsR2(i,:)=mboot';
    catch
        fprintf('bootstrap error with sample %0.0f;likely no. bins with counts <2\n',i);
        data.semR2(i)=nan;
    end

end

    if length(data.rep1ind)==length(data.rep2ind)

        combinedcounts=data.origbincounts(:,data.rep1ind)+data.origbincounts(:,data.rep2ind);
        combinedsigma=[];
        combinedmu=[];
        for i=1:length(combinedcounts(:,1))
            [m,s]=normfit(edges,[],[],combinedcounts(i,:));
            combinedsigma(end+1)=s;
            combinedmu(end+1)=m;
        end

        data.combinedsigma=combinedsigma;
        data.combinedmu=combinedmu;
    end
end
