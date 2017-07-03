addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
%% load data that are deemed "good", ones that have enough coverage
clear
g=load('YFSI_gooddata.mat');
gooddata=g.gooddata;

%% go through all sequences and convert every letter to single elements in cell array
loop1seq=cell(length(gooddata.seqs),1);
loop2seq=loop1seq;
for i=1:length(gooddata.seqs)
    s=gooddata.loop1{i};
    loop1seq{i}=s;
    s=gooddata.loop2{i};
    loop2seq{i}=s;
end
%% pull out seqs with specified loop lengths
l1len=[15:25 30 40 50 60];
l2len=7;

endslen=[2 3];
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
l1l2seqs=cell(length(l1l2),1);


for i=1:length(l1l2)
    l1=loop1seq{l1l2(i)};
    l1start=l1(1:endslen(1));
    l1end=l1((end-endslen(2)+1):end);
    s=sprintf('%s%s%s',l1start,l1end,loop2seq{l1l2(i)});
    l1l2seqs{i}=s;
end

outcsvcell=cell(length(l1l2),totallooplength+1);
for i=1:length(l1l2)
    for j=1:totallooplength
        s=l1l2seqs{i};
        outcsvcell{i,j}=s(j);
    end
    outcsvcell{i,totallooplength+1}=sprintf('%0.4f',gooddata.mus(l1l2(i)));
end

csvtitle=strcat(sprintf('l1+%0.0f,',1:endslen(1)),...
                sprintf('l1-%0.0f,',(endslen(2)):-1:1),...
                sprintf('l2+%0.0f,',1:l2len),sprintf('mu\n'));

fid = fopen('end3-3_N7.csv','w');
fprintf(fid,'%s\n',csvtitle);
for i=1:length(l1l2)
    s=strcat(sprintf('%s,',outcsvcell{i,1:end-1}),outcsvcell(i,end));
    fprintf(fid,'%s\n',s{1});
end
fclose(fid);

%% pull out seqs with specified loop lengths
l1len=7;
l2len=[15:25 30 40 50 60];

endslen=[3 3];
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
% mu25th=gooddata.mus<prctile|(gooddata.mus,25);
% l1l2=find(l1ind.*l2ind.*mu25th);

l1l2=find(l1ind.*l2ind);
l1l2seqs=cell(length(l1l2),1);


for i=1:length(l1l2)
    l2=loop2seq{l1l2(i)};
    l2start=l2(1:endslen(1));
    l2end=l2((end-endslen(2)+1):end);
    s=sprintf('%s%s%s',loop1seq{l1l2(i)},l2start,l2end);
    l1l2seqs{i}=s;
end

outcsvcell=cell(length(l1l2),totallooplength+1);
for i=1:length(l1l2)
    for j=1:totallooplength
        s=l1l2seqs{i};
        outcsvcell{i,j}=s(j);
    end
    outcsvcell{i,totallooplength+1}=sprintf('%0.4f',gooddata.mus(l1l2(i)));
end

csvtitle=strcat(sprintf('l1+%0.0f,',1:l1len),...
                sprintf('l2+%0.0f,',1:endslen(1)),...
                sprintf('l2-%0.0f,',(endslen(2)):-1:1),...
                sprintf('mu\n'));

fid = fopen('N7_end3-3.csv','w');
fprintf(fid,'%s\n',csvtitle);
for i=1:length(l1l2)
    s=strcat(sprintf('%s,',outcsvcell{i,1:end-1}),outcsvcell(i,end));
    fprintf(fid,'%s\n',s{1});
end
fclose(fid);
