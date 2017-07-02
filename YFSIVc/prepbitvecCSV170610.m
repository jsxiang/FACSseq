addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/FACSseq/
addpath ~/Documents/MATLAB/Standard/

%% load data that are deemed "good", ones that have enough coverage
% clear
% g=load('rbz.mat');
% DNA={'A','T','C','G'};
% 
% 
% data=findSEMbyBootstrp(g.rbz,100);
% %%
% save('rbz_bootstrp.mat','data');

g=load('rbz_bootstrp.mat');

data=g.data;
%%
motifs={'(GCTGTCACCGGA)([A|T|C|G]+)(TCCGGTCTGATGAGTCC)([A|T|C|G]+)(GGACGAAACAGC)',...
       '(GCTGTCACCGGA)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G]+)(TCCGGTCTGATGAGTCC)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G])(GGACGAAACAGC)',...
       '(GCTGTCACCGGA)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G])(TCCGGTCTGATGAGTCC)([A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G]+)(GGACGAAACAGC)'};
motifname={'sTRSV','smallLoop1','smallLoop2'};

sTRSVdata=findMotif(data,motifs,motifname);

%%
sTRSVdata(2).semC=data.semC(sTRSVdata(2).hasmotif);
sTRSVdata(3).semC=data.semC(sTRSVdata(3).hasmotif);
sTRSVdata(1).semC=data.semC(sTRSVdata(1).hasmotif);

%%
i=1
gooddata=sTRSVdata(i);
gooddata.standmus=Standard(gooddata.VYBmus,2,mean(gooddata.VYBmus),std(gooddata.VYBmus));
gooddata.standmus=mean(gooddata.standmus,2);

setfig('distribution');clf
h=histogram(gooddata.standmus)
h.NumBins=50

setfig('weights');clf

histogram(gooddata.semC)
set(gca,'XScale','log')

setfig('distribution vs weights');clf
plot(gooddata.standmus,1./(gooddata.semC.^2),'.','MarkerSize',1)
set(gca,'YScale','log')

%%

setfig('count vs weights');clf
plot(sum(gooddata.scaledbincounts,2),1./(gooddata.semC.^2),'.','MarkerSize',1)
% set(gca,'YScale','log')
set(gca,'XScale','log')

%%
ri=randperm(length(gooddata.seqs),length(gooddata.seqs));
gooddatarand.seqs=gooddata.seqs(ri);
gooddatarand.VYBmus=gooddata.standmus(ri);
gooddatarand.bincounts=gooddata.scaledbincounts(ri,:);
gooddatarand.weights=1./(gooddata.semC.^2);

traindata=struct;
traindata.seqs=gooddatarand.seqs(1:round(0.9*length(gooddatarand.seqs)));
traindata.VYBmus=gooddatarand.VYBmus(1:round(0.9*length(gooddatarand.seqs)));
traindata.bincounts=gooddatarand.bincounts(1:round(0.9*length(gooddatarand.seqs)),:);
traindata.weights=gooddatarand.weights(1:round(0.9*length(gooddatarand.seqs)));

% validdata=struct;
% validdata.seqs=gooddatarand.seqs((round(0.8*length(gooddatarand.seqs))+1):round(0.9*length(gooddatarand.seqs)));
% validdata.VYBmus=gooddatarand.VYBmus((round(0.8*length(gooddatarand.seqs))+1):round(0.9*length(gooddatarand.seqs)));
% validdata.bincounts=gooddatarand.bincounts((round(0.8*length(gooddatarand.seqs))+1):round(0.9*length(gooddatarand.seqs)),:);

testdata=struct;
testdata.seqs=gooddatarand.seqs((round(0.9*length(gooddatarand.seqs))+1):end);
testdata.VYBmus=gooddatarand.VYBmus((round(0.9*length(gooddatarand.seqs))+1):end);
testdata.bincounts=gooddatarand.bincounts((round(0.9*length(gooddatarand.seqs))+1):end,:);
testdata.weights=gooddatarand.weights((round(0.9*length(gooddatarand.seqs))+1):end);

% data=[validdata testdata traindata];
data=[testdata traindata];
outputfilename={'sTRSV_valid.mat','sTRSV_train.mat'};
%%
for k=1:length(data)
% plot enrichment versus mu
setfig('enrichment mu correlation');clf
plot(data(k).VYBmus,log2(sum(data(k).bincounts,2)),'.')
xlabel('\mu')
ylabel('log2(coverage)')


counts=log2(sum(data(k).bincounts,2));
muthresh=prctile(data(k).VYBmus,18);
mugood=data(k).VYBmus<muthresh;
% enrichedcounts=counts;
enrichedcounts=ones(length(data(k).seqs),1);

% enrichedcounts(mugood)=(enrichedcounts(mugood)-3); % arbitrarily enrich seqs according to their mu

% hold on
% plot(data(k).VYBmus,ceil(enrichedcounts),'.')
% hold off
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

% for i=1:length(data(k).seqs)
%     ec=ceil(enrichedcounts(i));
%     if data(k).VYBmus(i)<muthresh
%     expandedseqs=repmat(data(k).seqs(i),ec,1);
%     expandedmus=repmat(data(k).VYBmus(i),ec,1);
%     if length(enrichedseqs)>1
%         enrichedseqs={enrichedseqs{:},expandedseqs{:}};
%         enrichedmus=[enrichedmus(:);expandedmus];
%     else
%         enrichedseqs=expandedseqs;
%         enrichedmus=expandedmus;
%     end
%     else
%         enrichedseqs{end+1}=data(k).seqs{i};
%         enrichedmus(end+1)=data(k).VYBmus(i);
%     end
% end


% 
% for i=1:length(data(k).seqs)
%     ec=ceil(enrichedcounts(i));
%     if data(k).VYBmus(i)<muthresh
%     expandedmus=repmat(data(k).VYBmus(i),ec,1);
%     if length(enrichedseqs)>1
%         enrichedmus=[enrichedmus(:);expandedmus];
%     else
%         enrichedmus=expandedmus;
%     end
%     else
%         enrichedmus(end+1)=data(k).VYBmus(i);
%     end
% end



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




% %%
% addpath ~/Documents/MATLAB/Standard/
% enriched.standmus=Standard(enriched.mus,2,mean(enriched.mus),std(enriched.mus));
% 
% setfig('hist standmus enriched');clf
% histogram(enriched.standmus)

tr=struct;
tr.trainxdata=bitmat1;
tr.traindata=enriched.mus;

save(outputfilename{k},'tr');


end

%% divide sequences cleaver or non-cleaver down the middle
mus=g1.gooddata.mus;
switches1=mus<prctile(mus,50);

% mus=mean(g2.combined.cmus,2);
% switches2=mus<prctile(mus,20);
% switches=[switches1 switches2'];
% 
% % concatenate two structures
% gooddata.seqs=[g1.gooddata.seqs g2.combined.seqs];
% gooddata.loop1={g1.gooddata.loop1{:} g2.combined.loop1{:}};
% gooddata.loop2={g1.gooddata.loop2{:} g2.combined.loop2{:}};
% gooddata.loop1len=[g1.gooddata.loop1len g2.combined.loop1len];
% gooddata.loop2len=[g1.gooddata.loop2len g2.combined.loop2len];

switches=switches1;
gooddata=g1.gooddata;
%% go through all sequences and convert to matrix of 1 2 3 4
loop1seq=cell(length(gooddata.seqs),1);
loop2seq=loop1seq;
for i=1:length(gooddata.seqs)
    s=gooddata.loop1{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop1seq{i}=s-'0';
    
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


%% pull out seqs with specified loop lengths
l1len=8;
l2len=[15:25 30 40 50 60];

endslen=[7 7];
totallooplength=l1len+sum(endslen);
l1ind=zeros(1,length(gooddata.seqs));
for i=1:length(l1len)
    l1ind=l1ind| gooddata.loop1len==l1len(i);
end

l2ind=zeros(1,length(gooddata.seqs));
for i=1:length(l2len)
    l2ind=l2ind| gooddata.loop2len==l2len(i);
end


% let's apply an enrichment to only look at mu's less than the 25th
% percentile
% mu25th=gooddata.mus<prctile(gooddata.mus,25);
% l1l2=find(l1ind.*l2ind.*mu25th);
l1l2all=find(l1ind.*l2ind);
l1l2=find(l1ind.*l2ind);
l1l2seqs=zeros(length(l1l2),totallooplength);
for i=1:length(l1l2)
    l2=gooddata.loop2seqnum{l1l2(i)};
    l2start=l2(1:endslen(1));
    l2end=l2((end-endslen(2)+1):end);
    
    l1l2seqs(i,(l1len+1):(l1len+endslen(1)))=l2start;
    l1l2seqs(i,(l1len+endslen(1)+1):(l1len+endslen(1)+endslen(2)))=l2end;
   
    l1l2seqs(i,1:l1len)=gooddata.loop1seqnum{l1l2(i)};
end
l1l2seqsall=zeros(length(l1l2all),totallooplength);
for i=1:length(l1l2all)
    l2=gooddata.loop2seqnum{l1l2all(i)};
    l2start=l2(1:endslen(1));
    l2end=l2((end-endslen(2)+1):end);
    
    l1l2seqsall(i,(l1len+1):(l1len+endslen(1)))=l2start;
    l1l2seqsall(i,(l1len+endslen(1)+1):(l1len+endslen(1)+endslen(2)))=l2end;
   
    l1l2seqsall(i,1:l1len)=gooddata.loop1seqnum{l1l2all(i)};
end


%% Which positions to form 2-mer 
% try all possible 2-mers first
nmer=2;
twomerpos=zeros(1,nmer);
k=1;
for i=1:totallooplength
    for j=1:totallooplength
        if j<i
        twomerpos(k,1)=i;
        twomerpos(k,2)=j;
        k=k+1;
        end
    end
end

        


%% Incorporate 2-mer based upon loop-interactions


twomer=zeros(length(l1l2),length(twomerpos(:,1)));
feat2mer={};
for m=1:length(twomerpos(:,1))
for md=1:4
for nd=1:4
feat2mer{end+1}=sprintf('%0.0f-%0.0f:%s-%s',twomerpos(m,1),twomerpos(m,2),DNA{md},DNA{nd});
end
end
end

for i=1:length(l1l2)
    for m=1:length(twomerpos(:,1))
        for md=1:4
            for nd=1:4
                if l1l2seqs(i,twomerpos(m,1))==md && l1l2seqs(i,twomerpos(m,2))==nd
                    twomer(i,16*(m-1)+4*(md-1)+nd)=1;
                end
            end
        end
    end
end

setfig('two mer');clf
bar(sum(twomer))
%%
feats={};
for i=1:l1len
    for j=1:4
    feats{end+1}=sprintf('l1+%0.0f%s',i,DNA{j});
    end
end

for i=1:endslen(1)
    for j=1:4
    feats{end+1}=sprintf('l2+%0.0f%s',i,DNA{j});
    end
end

for i=(endslen(2)):-1:1
    for j=1:4
    feats{end+1}=sprintf('l2-%0.0f%s',i,DNA{j});
    end
end


feats={feats{:} feat2mer{:}};
feats{end+1}='mu';

% bitmat0=[bitmat twomer];
% 
% fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat_largeII.csv','w');
% csvtitle=sprintf('%s,',feats{:});
% fprintf(fid,'%s\n',csvtitle);
% 
% for i=1:length(l1l2)
%     s=sprintf('%0.0f,',bitmat0(i,1:end));
%     ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
%     fprintf(fid,'%s\n',ss);
% end
% fclose(fid);
% 
% fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat_largeII_2meronly.csv','w');
% csvtitle=sprintf('%s,',feats{end-length(feat2mer):end});
% fprintf(fid,'%s\n',csvtitle);
% 
% for i=1:length(l1l2)
%     s=sprintf('%0.0f,',twomer(i,:));
%     ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
%     fprintf(fid,'%s\n',ss);
% end
% fclose(fid);twomer2
%%
%% convert to bit 
bitmat2=zeros(length(l1l2), totallooplength,4);

for i=1:length(l1l2)
    for j=1:4
        for k=1:totallooplength
            bitmat2(i,k,j)=l1l2seqs(i,k)==j;
        end
    end
end

twomer2=zeros(length(l1l2),length(twomer(1,:))/16,4*4);
for i=1:length(l1l2)
    for j=1:16
        for k=1:(length(twomer(1,:))/16)
            twomer2(i,k,j)=twomer(i,(k-1)*16+j);
        end
    end
end

%% Train, valid, test
ri=randperm(length(twomer2),length(twomer2));
bitmatrand=twomer2(ri,:,:);

switchesrand=switches(l1l2(ri))';
tr=struct;
tr.trainxdata=bitmatrand(1:round(0.8*length(bitmatrand(:,1,1))),:,:);
tr.traindata=switchesrand(1:round(0.8*length(bitmatrand(:,1,1))));
save('/Volumes/smolke-lab$/Joy/YFSI/train2mer.mat','tr');

v=struct;
v.validxdata=bitmatrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))),:,:);
v.validdata=switchesrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))));
save('/Volumes/smolke-lab$/Joy/YFSI/valid2mer.mat','v');

tt=struct;
tt.testxdata=bitmatrand((round(0.9*length(bitmatrand(:,1,1)))+1):end,:,:);
tt.testdata=switchesrand((round(0.9*length(bitmatrand(:,1,1)))+1):end);
save('/Volumes/smolke-lab$/Joy/YFSI/test2mer.mat','tt');






%% END




