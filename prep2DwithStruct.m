addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/Standard/

% addpath ~/Documents/MATLAB/BREWER/
clear
%% load data that are deemed "good", ones that have enough coverage
g1=load('YFSI_gooddataWstruct.mat');
% g1=g1.gooddata;
load('MyColormaps','mycmap')
DNA={'A','T','C','G'};

validm=0.2156;
validc=-1.6998;

g1.gooddata.VYBmus=validm*g1.gooddata.mus+validc;
g1.gooddata.standmus=Standard(g1.gooddata.VYBmus,2,mean(g1.gooddata.VYBmus),std(g1.gooddata.VYBmus));
% filter for sequences that fold into ribozyme structure at all
structfilter = (cellfun('length',g1.gooddata.loop1struct)>0) & (cellfun('length',g1.gooddata.loop2struct)>0);


%% re-extract loop1struct and loop2struct
loop1structcell={};
loop2structcell={};
for i=1:length(g1.gooddata.ensembstruct)
    loop1structcell{end+1}=g1.gooddata.ensembstruct{i}(13:(12+g1.gooddata.loop1len(i)));
    loop2structcell{end+1}=g1.gooddata.ensembstruct{i}((end-11-g1.gooddata.loop2len(i)):(end-12));
end

g1.gooddata.loop1struct=loop1structcell;
g1.gooddata.loop2struct=loop2structcell;

%% save structure


%%
filter=(g1.gooddata.VYBmus<0.75) & structfilter;
gooddata.seqs=g1.gooddata.seqs(filter);
gooddata.VYBmus=g1.gooddata.VYBmus(filter);

% gooddata.VYBmus=(gooddata.VYBmus-min(gooddata.VYBmus))./(max(gooddata.VYBmus)-min(gooddata.VYBmus));
% gooddata.VYBmus=1-gooddata.VYBmus;
gooddata.loop1=g1.gooddata.loop1(filter);
gooddata.loop2=g1.gooddata.loop2(filter);
gooddata.loop1len=g1.gooddata.loop1len(filter);
gooddata.loop2len=g1.gooddata.loop2len(filter);
gooddata.loop1struct=g1.gooddata.loop1struct(filter);
gooddata.loop2struct=g1.gooddata.loop2struct(filter);
gooddata.ensembstruct=g1.gooddata.ensembstruct(filter);
gooddata.standmus=g1.gooddata.standmus(filter);

%% go through all sequences and convert to matrix of 1 2 3 4
loop1seq=cell(length(gooddata.seqs),1);
loop2seq=loop1seq;
loop1struct=loop1seq;
loop2struct=loop1seq;
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
    
    st=gooddata.loop1struct{i};
    st=regexprep(st,'(','1');
    st=regexprep(st,'\.','2');
    st=regexprep(st,')','3');
    loop1struct{i}=st-'0';
    
    st=gooddata.loop2struct{i};
    st=regexprep(st,'(','1');
    st=regexprep(st,'\.','2');
    st=regexprep(st,')','3');
    loop2struct{i}=st-'0';
end

gooddata.loop1seqnum=loop1seq;
gooddata.loop2seqnum=loop2seq;
gooddata.loop1structnum=loop1struct;
gooddata.loop2structnum=loop2struct;

%% pull out seqs with specified loop lengths
l2len=[15:25 30 40 50 60];
endslen=[7 7];
maxlen=max([l2len 4:8]);
l1l2seqs=[];
l1l2=[];
loop1seqs=[];
% loop2seqs=[];

muth=gooddata.VYBmus<prctile(gooddata.VYBmus,100); % filter
for k=4:8
    l1len=k;

    [seqs,index]=findSeqs(gooddata,l1len,l2len,endslen,muth);
%     loop1=seqs(:,1:k);
%     if length(loop1(1,:))<maxlen
%         loop1=[loop1 zeros(length(loop1(:,1)),maxlen-length(loop1(1,:)))];
%     end
%    
%     loop1seqs=[loop1seqs;loop1];
%     loop2=seqs(:,(k+1):end);
%     loop2seqs=[loop2seqs;loop2];
    l1l2=[l1l2 index];
end

loop2lenALL=gooddata.loop2len(l1l2);

loop2seqs=gooddata.loop2seqnum(l1l2);
halflen=floor(loop2lenALL/2);
loop2left=zeros(length(l1l2),max(halflen));
loop2right=zeros(length(l1l2),(max(halflen)));
for i=1:length(l1l2)
        loop2left(i,(end-length((halflen(i)+1):loop2lenALL(i))+1):end)=loop2seqs{i}((halflen(i)+1):loop2lenALL(i));
        loop2right(i,1:length(1:halflen(i)))=loop2seqs{i}(1:halflen(i));

end

loop2endjoin=[loop2left loop2right];


loop1seqs=gooddata.loop1seqnum(l1l2);
loop1lenALL=gooddata.loop1len(l1l2);
halflen=floor(loop1lenALL/2);
loop1left=zeros(length(l1l2),max(halflen));
loop1right=zeros(length(l1l2),max(halflen));
for i=1:length(l1l2)
    loop1left(i,(end-length((halflen(i)+1):loop1lenALL(i))+1):end)=loop1seqs{i}((halflen(i)+1):loop1lenALL(i));
    loop1right(i,1:length(1:halflen(i)))=loop1seqs{i}(1:halflen(i));
end

loop1endjoin=[loop1left loop1right];
% center and invert
lendiff=abs(length(loop2endjoin(1,:))-length(loop1endjoin(1,:)));
loop1endjoin=[zeros(length(l1l2),lendiff/2) loop1endjoin zeros(length(l1l2),lendiff/2)];
loop1endjoin=loop1endjoin(:,end:-1:1);

%% pull out structures with specified loop lengths
muth=gooddata.VYBmus<prctile(gooddata.VYBmus,100); % filter
loop1structures=[];

loop2lenALL=gooddata.loop2len(l1l2);

loop2structures=gooddata.loop2structnum(l1l2);
halflen=floor(loop2lenALL/2);
loop2left=zeros(length(l1l2),max(halflen));
loop2right=zeros(length(l1l2),(max(halflen)));
for i=1:length(l1l2)
        loop2left(i,(end-length((halflen(i)+1):loop2lenALL(i))+1):end)=loop2structures{i}((halflen(i)+1):loop2lenALL(i));
        loop2right(i,1:length(1:halflen(i)))=loop2structures{i}(1:halflen(i));

end

loop2endjoin_structures=[loop2left loop2right];


loop1structures=gooddata.loop1structnum(l1l2);
loop1lenALL=gooddata.loop1len(l1l2);
halflen=floor(loop1lenALL/2);
loop1left=zeros(length(l1l2),max(halflen));
loop1right=zeros(length(l1l2),max(halflen));
for i=1:length(l1l2)
    loop1left(i,(end-length((halflen(i)+1):loop1lenALL(i))+1):end)=loop1structures{i}((halflen(i)+1):loop1lenALL(i));
    loop1right(i,1:length(1:halflen(i)))=loop1structures{i}(1:halflen(i));
end

loop1endjoin_structures=[loop1left loop1right];
% center and invert
lendiff=abs(length(loop2endjoin_structures(1,:))-length(loop1endjoin_structures(1,:)));
loop1endjoin_structures=[zeros(length(l1l2),lendiff/2) loop1endjoin_structures zeros(length(l1l2),lendiff/2)];
loop1endjoin_structures=loop1endjoin_structures(:,end:-1:1);




%% convert ensembstruct to one hot encoding

%% convert to bit 8 pixels by 1 channels
l1l2seqs(:,:,1)=loop1endjoin;
l1l2seqs(:,:,2)=loop2endjoin;

l1l2structs(:,:,1)=loop1endjoin_structures;
l1l2structs(:,:,2)=loop2endjoin_structures;
bitmat2=zeros(length(l1l2),maxlen,2,8);

for i=1:length(l1l2)
    for j=1:maxlen
        for q=1:2
            for p=1:4
                bitmat2(i,j,1,4*(q-1)+p)=l1l2seqs(i,j,q)==p;
            end
            for r=1:3
                bitmat2(i,j,2,4*(q-1)+r)=l1l2structs(i,j,q)==r;
            end
        end

    end
end




%% Train, valid, test
switches=gooddata.standmus(l1l2);

ri=randperm(length(bitmat2),length(bitmat2));
bitmatrand=bitmat2(ri,:,:,:);
switchesrand=switches(ri)';
tr=struct;

bitmatlen=length(bitmatrand(:,1,1,1));

tr.trainxdata=bitmatrand(1:round(0.95*bitmatlen),:,:,:);
tr.traindata=switchesrand(1:round(0.95*bitmatlen));
save('~/Documents/rbznn/train4DwStruct.mat','tr');

v=struct;
v.validxdata=bitmatrand((round(0.95*bitmatlen)+1):end,:,:,:);
v.validdata=switchesrand((round(0.95*bitmatlen)+1):end);
save('~/Documents/rbznn/valid4DwStruct.mat','v');

% tt=struct;
% tt.testxdata=bitmatrand((round(0.9*bitmatlen)+1):end,:,:);
% tt.testdata=switchesrand((round(0.9*bitmatlen)+1):end);
% save('~/Documents/rbznn/test3D.mat','tt');
% 









%% END














