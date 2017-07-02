addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
%% load data that are deemed "good", ones that have enough coverage
% clear
g1=load('YFSI_gooddata.mat');
g2=load('~/Documents/YFS/YFSIII/YFSII_gooddata.mat');
gooddata=g1.gooddata;
DNA={'A','T','C','G'};

%% go through all sequences and convert to matrix of 1 2 3 4
% note: full length sequences are considered
seq=cell(length(gooddata.seqs),1);
% loop2seq=loop1seq;
for i=1:length(gooddata.seqs)
    s=gooddata.seqs{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    seq{i}=s-'0';
%     
%     s=gooddata.loop2{i};
%     s=regexprep(s,'A','1');
%     s=regexprep(s,'T','2');
%     s=regexprep(s,'C','3');
%     s=regexprep(s,'G','4');
%     % this is a super cool way of converting a string of numbers to a matrix
%     loop2seq{i}=s-'0';
%     
end

gooddata.seqnum=seq;


%% pull out seqs with specified loop lengths
% l1len=[15:25 30 40 50 60];
% l2len=8;
l1len=[4 5 6 7 8];
l1len=7;
l2len=[15:max(gooddata.loop2len)];

% endslen=[15 15];
% totallooplength=sum(endslen)+l2len;
l1ind=zeros(1,length(gooddata.seqs));   
for i=1:length(l1len)
    l1ind=l1ind|gooddata.loop1len==l1len(i); % This "OR" is serving as an addition
end

l2ind=zeros(1,length(gooddata.seqs));
for i=1:length(l2len)
    l2ind=l2ind|gooddata.loop2len==l2len(i);
end
    
% let's apply an enrichment to only look at mu's less than the 25th
% percentile
% mu25th=gooddata.mus<prctile(gooddata.mus,25);
% l1l2=find(l1ind.*l2ind.*mu25th);

l1l2=find(l1ind.*l2ind);
rbzlen=41; % backbone ribozyme length
totlen=max(l1len)+max(l2len)+rbzlen;
l1l2seqs=zeros(length(l1l2),totlen);


for i=1:length(l1l2)
    l1l2seqs(i,1:length(gooddata.seqnum{l1l2(i),:}))=gooddata.seqnum{l1l2(i),:};
    % data is already zero padded
end
%% convert to bit 
bitmat1=zeros(length(l1l2), totlen,4);

for i=1:length(l1l2)
    for j=1:4
        for k=1:totlen
            bitmat1(i,k,j)=l1l2seqs(i,k)==j;
        end
    end
end


%% divide sequences cleaver or non-cleaver down the middle
mus=gooddata.mus(l1l2);
switches1=mus<prctile(mus,15);


%%
g2=load('~/Documents/YFS/YFSIII/YFSII_gooddata.mat');
gooddata=g2.combined;
DNA={'A','T','C','G'};

%% go through all sequences and convert to matrix of 1 2 3 4
% note: full length sequences are considered
seq=cell(length(gooddata.seqs),1);
% loop2seq=loop1seq;
for i=1:length(gooddata.seqs)
    s=gooddata.seqs{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    seq{i}=s-'0';
%     
%     s=gooddata.loop2{i};
%     s=regexprep(s,'A','1');
%     s=regexprep(s,'T','2');
%     s=regexprep(s,'C','3');
%     s=regexprep(s,'G','4');
%     % this is a super cool way of converting a string of numbers to a matrix
%     loop2seq{i}=s-'0';
%     
end

gooddata.seqnum=seq;


%% pull out seqs with specified loop lengths
% l1len=[15:25 30 40 50 60];
% l2len=8;
l1len=[4 5 6 7 8];
l1len=7;
l2len=[15:max([gooddata.loop2len l2len])];

% endslen=[15 15];
% totallooplength=sum(endslen)+l2len;
l1ind=zeros(1,length(gooddata.seqs));   
for i=1:length(l1len)
    l1ind=l1ind|gooddata.loop1len==l1len(i); % This "OR" is serving as an addition
end

l2ind=zeros(1,length(gooddata.seqs));
for i=1:length(l2len)
    l2ind=l2ind|gooddata.loop2len==l2len(i);
end
    
% let's apply an enrichment to only look at mu's less than the 25th
% percentile
% mu25th=gooddata.mus<prctile(gooddata.mus,25);
% l1l2=find(l1ind.*l2ind.*mu25th);

l1l2=find(l1ind.*l2ind);
rbzlen=41; % backbone ribozyme length
totlen=max(l1len)+max(l2len)+rbzlen;
l1l2seqs=zeros(length(l1l2),totlen);


for i=1:length(l1l2)
    l1l2seqs(i,1:length(gooddata.seqnum{l1l2(i),:}))=gooddata.seqnum{l1l2(i),:};
    % data is already zero padded
end
%% convert to bit 
bitmat2=zeros(length(l1l2), totlen,4);

for i=1:length(l1l2)
    for j=1:4
        for k=1:totlen
            bitmat2(i,k,j)=l1l2seqs(i,k)==j;
        end
    end
end

%% divide sequences cleaver or non-cleaver down the middle
mus=mean(gooddata.cmus(l1l2,:),2);
switches2=mus<prctile(mus,20);
bitmat2=[];
bitmat=[bitmat1;bitmat2];
switches=[switches1 switches2'];
switches=switches1;

%% Train, valid, test
ri=randperm(length(bitmat),length(bitmat));
bitmatrand=bitmat(ri,:,:);
switchesrand=switches(ri)';
tr=struct;
tr.trainxdata=bitmatrand(1:round(0.8*length(bitmatrand(:,1,1))),:,:);
tr.traindata=switchesrand(1:round(0.8*length(bitmatrand(:,1,1))));
save('~/Documents/CS273B/trainlstm.mat','tr');

v=struct;
v.validxdata=bitmatrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))),:,:);
v.validdata=switchesrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))));
save('~/Documents/CS273B/validlstm.mat','v');

tt=struct;
tt.testxdata=bitmatrand((round(0.9*length(bitmatrand(:,1,1)))+1):end,:,:);
tt.testdata=switchesrand((round(0.9*length(bitmatrand(:,1,1)))+1):end);
save('~/Documents/CS273B/testlstm.mat','tt');



%% Alternatively, include all sequences with l1len = 8 and l2len = all
% 
% %% pull out seqs with specified loop lengths
% % l1len=[15:25 30 40 50 60];
% % l2len=8;
% l1len=8;
% l2len=[15:25 30 40 50 60];
% 
% % endslen=[15 15];
% % totallooplength=sum(endslen)+l2len;
% l1ind=zeros(1,length(gooddata.seqs));   
% for i=1:length(l1len)
%     l1ind=l1ind|gooddata.loop1len==l1len(i); % This "OR" is serving as an addition
% end
% 
% l2ind=zeros(1,length(gooddata.seqs));
% for i=1:length(l2len)
%     l2ind=l2ind|gooddata.loop2len==l2len(i);
% end
%     
% % let's apply an enrichment to only look at mu's less than the 25th
% % percentile
% % mu25th=gooddata.mus<prctile(gooddata.mus,25);
% % l1l2=find(l1ind.*l2ind.*mu25th);
% 
% l1l2=find(l1ind.*l2ind);
% rbzlen=41; % backbone ribozyme length
% totlen=l1len+max(l2len)+rbzlen;
% l1l2seqs=zeros(length(l1l2),totlen);
% 
% 
% for i=1:length(l1l2)
%     l1l2seqs(i,1:length(gooddata.seqnum{l1l2(i),:}))=gooddata.seqnum{l1l2(i),:};
%     % data is already zero padded
% end
% %% divide sequences cleaver or non-cleaver down the middle
% setfig('mus');clf
% hist(gooddata.mus(l1l2),150)
% mus=gooddata.mus(l1l2);
% switches=mus<prctile(mus,20);
% 
% %% convert to bit 
% bitmat=zeros(length(l1l2), totlen,4);
% 
% for i=1:length(l1l2)
%     for j=1:4
%         for k=1:totlen
%             bitmat(i,k,j)=l1l2seqs(i,k)==j;
%         end
%     end
% end
% 
% 
% %% Train, valid, test
% ri=randperm(length(bitmat),length(bitmat));
% bitmatrand=bitmat(ri,:,:);
% switchesrand=switches(ri)';
% tr=struct;
% tr.trainxdata=bitmatrand(1:round(0.8*length(bitmatrand(:,1,1))),:,:);
% tr.traindata=switchesrand(1:round(0.8*length(bitmatrand(:,1,1))));
% save('~/Documents/CS273B/trainlstm.mat','tr');
% 
% v=struct;
% v.validxdata=bitmatrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))),:,:);
% v.validdata=switchesrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))));
% save('~/Documents/CS273B/validlstm.mat','v');
% 
% tt=struct;
% tt.testxdata=bitmatrand((round(0.9*length(bitmatrand(:,1,1)))+1):end,:,:);
% tt.testdata=switchesrand((round(0.9*length(bitmatrand(:,1,1)))+1):end);
% save('~/Documents/CS273B/testlstm.mat','tt');



%% END


