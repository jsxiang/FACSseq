addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
%% load data that are deemed "good", ones that have enough coverage
% clear
g=load('YFSI_gooddata.mat');
gooddata=g.gooddata;
DNA={'A','T','C','G'};

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
% l1len=[15:25 30 40 50 60];
% l2len=8;
l1len=[30];
l2len=7;

endslen=[3 3];
totallooplength=sum(endslen)+l2len;
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

l1l2=find(l1ind.*l2ind);
l1l2seqs=zeros(length(l1l2),totallooplength);
for i=1:length(l1l2)
    l1=gooddata.loop1seqnum{l1l2(i)};
    l1start=l1(1:endslen(1));
    l1end=l1((end-endslen(2)+1):end);
    
    l1l2seqs(i,1:endslen(1))=l1start;
    l1l2seqs(i,(endslen(1)+1):(endslen(1)+endslen(2)))=l1end;
   
    l1l2seqs(i,(sum(endslen)+1):(sum(endslen)+l2len))=gooddata.loop2seqnum{l1l2(i)};
end

%% convert to bit 
bitmat=zeros(length(l1l2),length(l1l2seqs(1,:))*4);

for i=1:4
    bitmat(:,i:4:end)=l1l2seqs==i;
end

%% Incorporate 2-mer based upon loop-interactions
% run analyzePSC first to get consensus (see bottom of that script)
consensus =[
     5     1     3     2
     6     1     3     2
    13     1     4     2
     6     5     3     3
    13     5     4     4
    12     6     1     3
    13     6     4     2];

twomer=zeros(length(l1l2),length(consensus(:,1)));
feat2mer={};
for m=1:length(consensus(:,1))
for md=1:4
for nd=1:4
feat2mer{end+1}=sprintf('%0.0f-%0.0f:%s-%s',consensus(m,1),consensus(m,2),DNA{md},DNA{nd});
end
end
end

for i=1:length(l1l2)
    for m=1:length(consensus(:,1))
        for md=1:4
            for nd=1:4
                if l1l2seqs(i,consensus(m,1))==md && l1l2seqs(i,consensus(m,2))==nd
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
for i=1:endslen(1)
    for j=1:4
    feats{end+1}=sprintf('l1+%0.0f%s',i,DNA{j});
    end
end

for i=(endslen(2)):-1:1
    for j=1:4
    feats{end+1}=sprintf('l1-%0.0f%s',i,DNA{j});
    end
end

for i=1:l2len
    for j=1:4
    feats{end+1}=sprintf('l2+%0.0f%s',i,DNA{j});
    end
end
feats={feats{:} feat2mer{:}};
feats{end+1}='mu';

bitmat=[bitmat twomer];

fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat.csv','w');
csvtitle=sprintf('%s,',feats{:});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',bitmat(i,1:end));
    ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end
fclose(fid);

fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat_2meronly.csv','w');
csvtitle=sprintf('%s,',feats{end-length(feat2mer):end});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',twomer(i,:));
    ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end 
fclose(fid);




%% pull out seqs with specified loop lengths
l1len=8;
l2len=[15:25 30 40 50 60];

endslen=[4 4];
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
mu25th=gooddata.mus<prctile(gooddata.mus,25);
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


%% Incorporate 2-mer based upon loop-interactions
% run analyzePSC first to get consensus (see bottom of that script)


twomer=zeros(length(l1l2),length(consensus(:,1)));
feat2mer={};
for m=1:length(consensus(:,1))
for md=1:4
for nd=1:4
feat2mer{end+1}=sprintf('%0.0f-%0.0f:%s-%s',consensus(m,1),consensus(m,2),DNA{md},DNA{nd});
end
end
end

for i=1:length(l1l2)
    for m=1:length(consensus(:,1))
        for md=1:4
            for nd=1:4
                if l1l2seqs(i,consensus(m,1))==md && l1l2seqs(i,consensus(m,2))==nd
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

bitmat=[bitmat twomer];

fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat_largeII.csv','w');
csvtitle=sprintf('%s,',feats{:});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',bitmat(i,1:end));
    ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end
fclose(fid);

fid= fopen('/Volumes/smolke-lab$/Joy/YFSI/bitmat_largeII_2meronly.csv','w');
csvtitle=sprintf('%s,',feats{end-length(feat2mer):end});
fprintf(fid,'%s\n',csvtitle);

for i=1:length(l1l2)
    s=sprintf('%0.0f,',twomer(i,:));
    ss=sprintf('%s%0.4f',s,gooddata.mus(l1l2(i)));
    fprintf(fid,'%s\n',ss);
end
fclose(fid);


%% END


