addpath ../FACSseqlib/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACSseq/

addpath ~/Documents/MATLAB/FACS/
addpath ../FACSseqlib/FACSseq/
addpath ../FACSseqlib/Standard/

%% load data that are deemed "good", ones that have enough coverage
% clear
% g=load('rbz.mat');
% DNA={'A','T','C','G'};
% data=findSEMbyBootstrp(g.rbz,100);
% %%
% save('rbz_bootstrp.mat','data');

g=load('~/Documents/YFS/YFSIVc/rbz_bootstrp.mat');

data=g.data;

%%
motifs={'(^GCTGTCACCGGA)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G]+)(TCCGGTCTGATGAGTCC)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G]+)(GGACGAAACAGC)',...
       '(^GCTGTCACCGGA)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G]+)(TCCGGTCTGATGAGTCC)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G])(GGACGAAACAGC)',...
       '(^GCTGTCACCGGA)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G])(TCCGGTCTGATGAGTCC)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G]+)(GGACGAAACAGC)'};
motifname={'sTRSV','smallLoop1','smallLoop2'};

sTRSVdata=findMotif(data,motifs,motifname);


sTRSVdata(2).semC=data.semC(sTRSVdata(2).hasmotif);
sTRSVdata(3).semC=data.semC(sTRSVdata(3).hasmotif);
sTRSVdata(1).semC=data.semC(sTRSVdata(1).hasmotif);
sTRSVdata(1).semR1=data.semR1(sTRSVdata(1).hasmotif);
sTRSVdata(2).semR1=data.semR1(sTRSVdata(2).hasmotif);
sTRSVdata(1).semR2=data.semR2(sTRSVdata(1).hasmotif);
sTRSVdata(2).semR2=data.semR2(sTRSVdata(2).hasmotif);

%%
i=1
numseqs=50;
goodnum=find(sum(sTRSVdata(i).origbincounts,2)>numseqs);
gooddata.seqs=sTRSVdata(i).seqs(goodnum);

gooddata.VYBmus=sTRSVdata(i).VYBmus(goodnum,:);
gooddata.standmus=Standard(gooddata.VYBmus,2,mean(gooddata.VYBmus),std(gooddata.VYBmus));
gooddata.standmus=mean(gooddata.standmus,2);
gooddata.semC=sTRSVdata(i).semC(goodnum);
gooddata.semR1=sTRSVdata(i).semR1(goodnum);
gooddata.semR2=sTRSVdata(i).semR2(goodnum);
gooddata.origbincounts=sTRSVdata(i).origbincounts(goodnum,:);
%%
setfig('weights');clf

histogram(gooddata.semC)
set(gca,'XScale','log')

setfig('distribution vs weights');clf
plot(gooddata.standmus,1./(gooddata.semC.^2),'.','MarkerSize',1)
set(gca,'YScale','log')

setfig('count vs weights');clf
plot(sum(gooddata.origbincounts,2),1./(gooddata.semC.^2),'.','MarkerSize',1)
% set(gca,'YScale','log')
set(gca,'XScale','log')




%%
% assign bin using edges from histogram
% random sample without replacement, keep 750 per bin
setfig('distribution');clf
numbins=50;
gooddata.meanVYBmus=mean(gooddata.VYBmus,2);
h=histogram(gooddata.meanVYBmus,numbins);

%%
binthresh=7000;
evendata.seqs={};
evendata.VYBmus=[];
evendata.standmus=[];
evendata.weights=[];
evendata.weightsreps=[];

for i=1:(length(h.BinEdges)-1)
    bi=find(gooddata.meanVYBmus>(h.BinEdges(i)) & gooddata.meanVYBmus<(h.BinEdges(i+1)));
    if length(bi)<=binthresh
        binidx=bi;
    else
        binidx=datasample(bi,binthresh);

    end
    evendata.seqs={evendata.seqs{:},gooddata.seqs{binidx}};
    evendata.VYBmus=[evendata.VYBmus; gooddata.VYBmus(binidx,:)];
    evendata.weights=[evendata.weights; 1./(gooddata.semC(binidx).^2)];
    evendata.weightsreps=[evendata.weightsreps; 1./(gooddata.semR1(binidx).*gooddata.semR2(binidx))];

end
evendata.weightsreps(find(isinf(evendata.weightsreps)))=0;
evendata.weightsreps(find(isnan(evendata.weightsreps)))=0;

setfig('evenly distributed data');clf
dscatter(evendata.VYBmus(:,1),evendata.VYBmus(:,2),'bins',[50 50])

mdl = fitlm(evendata.VYBmus(:,1),evendata.VYBmus(:,2),'weights',evendata.weightsreps);
rsq=mdl.Rsquared.Adjusted;
text(-.5,2,sprintf('R^2 = %0.2f',rsq))


evendata.standmus=Standard(evendata.VYBmus,2,mean(evendata.VYBmus),std(evendata.VYBmus));
setfig('distribution');
hold on
histogram(mean(evendata.VYBmus,2),numbins)

histogram(evendata.standmus,numbins)

%%
ri=randperm(length(evendata.seqs),length(evendata.seqs));
gooddatarand.seqs=evendata.seqs(ri);
gooddatarand.VYBmus=evendata.standmus(ri,:);
gooddatarand.weights=evendata.weights;
%%
% save('gooddatarand.mat','gooddatarand');
%%
traindata=struct;
traindata.seqs=gooddatarand.seqs(1:round(0.9*length(gooddatarand.seqs)));
traindata.seqs={traindata.seqs{:} traindata.seqs{:}};
traindata.VYBmus=[gooddatarand.VYBmus(1:round(0.9*length(gooddatarand.seqs)),1);gooddatarand.VYBmus(1:round(0.9*length(gooddatarand.seqs)),2)];
% traindata.bincounts=gooddatarand.bincounts(1:round(0.9*length(gooddatarand.seqs)),:);
traindata.weights=gooddatarand.weights(1:round(0.9*length(gooddatarand.seqs)));

testdata=struct;
testdata.seqs=gooddatarand.seqs((round(0.9*length(gooddatarand.seqs))+1):end);
% testdata.seqs={testdata.seqs{:} testdata.seqs{:}};
testdata.VYBmus=mean(gooddatarand.VYBmus((round(0.9*length(gooddatarand.seqs))+1):end,:),2);
% testdata.bincounts=gooddatarand.bincounts((round(0.9*length(gooddatarand.seqs))+1):end,:);
testdata.weights=gooddatarand.weights((round(0.9*length(gooddatarand.seqs))+1):end);

% data=[validdata testdata traindata];
data=[testdata traindata];
outputfilename={'~/Documents/rbznn/sTRSVdupALL_valid.mat','~/Documents/rbznn/sTRSVdupALL_train.mat'};
% outputfilename={'~/Documents/rbznn/sTRSV_valid.mat','~/Documents/rbznn/sTRSV_train.mat'};
%%
for k=1:length(data)



muthresh=prctile(data(k).VYBmus,18);
mugood=data(k).VYBmus<muthresh;
enrichedcounts=ones(length(data(k).seqs),1);


enrichedseqs={};
enrichedmus=[];
enrichedweights=[];
for i=1:length(data(k).seqs)
enrichedseqs{end+1}=data(k).seqs{i};
enrichedmus(end+1)=data(k).VYBmus(i);

end

data(k).weights(isinf(log10(data(k).weights)))=0;

enrichedweights=log10(data(k).weights)./max(log10(data(k).weights));
enrichedweights(isinf(abs(enrichedweights))|isnan(enrichedweights))=0;

% go through all sequences and convert to matrix of 1 2 3 4
seq=cell(length(enrichedseqs),1);

for i=1:length(enrichedseqs)
    s=enrichedseqs{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    seq{i}=s-'0';
    

end


enriched.seqs=enrichedseqs;
enriched.mus=enrichedmus';
enriched.seqnum=seq;
enriched.weights=(enrichedweights);
%
totlen=113;
l1l2seqs=zeros(length(enriched.seqs),totlen);

for i=1:length(enriched.seqnum)
    l1l2seqs(i,1:length(enriched.seqnum{i}))=enriched.seqnum{i};
end


% convert to bit 
bitmat1=zeros(length(enriched.seqs), totlen,4);

for i=1:length(enriched.seqs)
    for j=1:4
        for l=1:totlen
            bitmat1(i,l,j)=(l1l2seqs(i,l)==j);
        end
    end
end

% consider weighting by 1/sigma
tr=struct;
tr.trainxdata=bitmat1;
tr.traindata=enriched.mus;

save(outputfilename{k},'tr');


end






%% END




