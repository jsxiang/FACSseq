addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%%

g=load('YFSI_gooddata.mat');
gooddata=g.gooddata;


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

%% for all seqs with specified loop lengths
% what is the distribution of mus
l1len=40;
l2len=8;
totallooplength=l1len+l2len;
l1l2=find(gooddata.loop1len==l1len&gooddata.loop2len==l2len);
l1l2seqs=zeros(length(l1l2),totallooplength);
for i=1:length(l1l2)
    l1l2seqs(i,1:l1len)=gooddata.loop1seqnum{l1l2(i)};
    l1l2seqs(i,(l1len+1):(l1len+l2len))=gooddata.loop2seqnum{l1l2(i)};
end



%% Analyses
findParameters(gooddata.mus,l1l2seqs);

% find mutual sequence contribution - pairwise mu
pairwise=findPairwiseMu(gooddata.mus(l1l2),l1l2seqs);

% find entropy
ent=findEnt(l1l2seqs);

% find mutual information
[mI,pJoint]=findMutualInformation(l1l2seqs);

%% END