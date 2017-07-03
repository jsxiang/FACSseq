clear
addpath ~/Documents/robot/Matlab-Utilities/

%% load in the sequences that occurred more than once
fid=fopen('../libreads/MFSI_libnoBC_MoreThanOnce.txt');
seqs=textscan(fid,'%s');
fclose(fid);

dat=load('MFSI.mat'); % this is loaded in as a struct

%%
totalcounts_perseq=sum(dat.MFSIdat,2);
totalcounts_perbarcode=sum(dat.MFSIdat);

dlmwrite('barcoded_count.txt',totalcounts_perbarcode')
%%
setfig('total counts');clf
hist(totalcounts_perseq,1000)
% set(gca,'YScale','log')
ylabel('Frequency')
xlabel('NGS read count')

%%
x=load('barcodeglobalcount.txt');
scalex=(x(1:40,1)./x(1:40,3))';
scalexmat=repmat(scalex,length(seqs{:}),1);
scaleddata=dat.MFSIdat(:,1:40).*scalexmat;
dat.seqs=seqs{:};
%%
dat.cond(1).scaled=scaleddata(:,1:10);
dat.cond(2).scaled=scaleddata(:,11:20);
dat.cond(3).scaled=scaleddata(:,21:30);
dat.cond(4).scaled=scaleddata(:,31:40);

cond1=dat.MFSIdat(:,1:10);
cond2=dat.MFSIdat(:,11:20);
cond3=dat.MFSIdat(:,21:30);
cond4=dat.MFSIdat(:,31:40);

%%
n=8;
goodcoveragerep1=(sum(cond1,2)>n) .* (sum(cond2,2)>n);
goodcoveragerep2=(sum(cond3,2)>n) .* (sum(cond4,2)>n);
goodcoverageALL=(sum(cond1,2)>n) .* (sum(cond2,2)>n) .* (sum(cond3,2)>n) .* (sum(cond4,2)>n);
goodcoverageCombined=(sum(cond1+cond3,2)>n) .* (sum(cond2+cond3,2)>n);
gooddata=scaleddata(find(goodcoverageALL),:);
goodseq=dat.seqs(find(goodcoverageALL));

%%
for i=1:length(goodseq)
    for j=1:4
        ind=(1+(j-1)*10):(j*10);
        goodcoverage(i).cond(j).data=gooddata(i,ind);
    end
    goodcoverage(i).seq=goodseq(i);
end

save('MFSI_goodcoverage.mat','goodcoverage')


gooddata=scaleddata(find(goodcoveragerep1),:);
goodseq=dat.seqs(find(goodcoveragerep1));

for i=1:length(goodseq)
    for j=1:4
        ind=(1+(j-1)*10):(j*10);
        goodcoverage(i).cond(j).data=gooddata(i,ind);
    end
    goodcoverage(i).seq=goodseq(i);
end
save('MFSI_goodcoveragerep1','goodcoverage')

gooddata=scaleddata(find(goodcoveragerep2),:);
goodseq=dat.seqs(find(goodcoveragerep2));

for i=1:length(goodseq)
    for j=1:4
        ind=(1+(j-1)*10):(j*10);
        goodcoverage(i).cond(j).data=gooddata(i,ind);
    end
    goodcoverage(i).seq=goodseq(i);
end
save('MFSI_goodcoveragerep2','goodcoverage')

gooddata=scaleddata(find(goodcoverageCombined),:);
goodseq=dat.seqs(find(goodcoverageCombined));

for i=1:length(goodseq)
    for j=1:4
        ind=(1+(j-1)*10):(j*10);
        goodcoverage(i).cond(j).data=gooddata(i,ind);
    end
    goodcoverage(i).seq=goodseq(i);
end
save('MFSI_goodcoverageCombined','goodcoverage')



%% REALLY good coverage
n=50;
goodcoveragerep1=(sum(cond1,2)>n) .* (sum(cond2,2)>n);
goodcoveragerep2=(sum(cond3,2)>n) .* (sum(cond4,2)>n);
goodcoverageALL=(sum(cond1,2)>n) .* (sum(cond2,2)>n) .* (sum(cond3,2)>n) .* (sum(cond4,2)>n);
goodcoverageCombined=(sum(cond1+cond3,2)>n) .* (sum(cond2+cond3,2)>n);
gooddata=scaleddata(find(goodcoverageALL),:);
goodseq=dat.seqs(find(goodcoverageALL));
for i=1:length(goodseq)
    for j=1:4
        ind=(1+(j-1)*10):(j*10);
        goodcoverage(i).cond(j).data=gooddata(i,ind);
    end
    goodcoverage(i).seq=goodseq(i);
end

save('MFSI_REALLYgoodcoverage.mat','goodcoverage')














