clear
addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/Standard/

% addpath ~/Documents/MATLAB/BREWER/
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
%%
structcell={};
for i=1:length(g1.gooddata.ensembstruct)
    structcell{end+1}=g1.gooddata.ensembstruct{i}(:);
end

%% 
f=fopen('goodseqs.txt','w');
fprintf(f,'%s\n',g1.gooddata.seqs{:});
fclose(f);

gooddata.standmus=Standard(gooddata.VYBmus,2,mean(gooddata.VYBmus),std(gooddata.VYBmus));

g=fopen('goodmus.txt','w');
fprintf(g,'%0.4f\n',g1.gooddata.standmus);
fclose(g);

h=fopen('goodstructs.txt','w');
for i=1:length(structcell)
    fprintf(h,'%s\n',structcell{1});
end
fclose(h);

%% save structure


%%
filter=(g1.gooddata.standmus<0.75) & structfilter;
gooddata.seqs=g1.gooddata.seqs(filter);
gooddata.VYBmus=g1.gooddata.VYBmus(filter);
gooddata.bincounts=g1.gooddata.bincounts(filter,:);

gooddata.ensembstruct=g1.gooddata.ensembstruct(filter);
gooddata.standmus=g1.gooddata.standmus(filter);
%%
ri=randperm(length(gooddata.seqs),length(gooddata.seqs));
gooddatarand.seqs=gooddata.seqs(ri);
gooddatarand.VYBmus=gooddata.standmus(ri);
gooddatarand.bincounts=gooddata.bincounts(ri,:);
gooddatarand.structs=gooddata.ensembstruct(ri);

%%
traindata=struct;
traindata.seqs=gooddatarand.seqs(1:round(0.95*length(gooddatarand.seqs)));
traindata.VYBmus=gooddatarand.VYBmus(1:round(0.95*length(gooddatarand.seqs)));
traindata.bincounts=gooddatarand.bincounts(1:round(0.95*length(gooddatarand.seqs)),:);
traindata.structs=gooddatarand.structs(1:round(0.95*length(gooddatarand.seqs)));

validdata=struct;
validdata.seqs=gooddatarand.seqs((round(0.95*length(gooddatarand.seqs))+1):end);
validdata.VYBmus=gooddatarand.VYBmus((round(0.95*length(gooddatarand.seqs))+1):end);
validdata.bincounts=gooddatarand.bincounts((round(0.95*length(gooddatarand.seqs))+1):end,:);
validdata.structs=gooddatarand.structs((round(0.95*length(gooddatarand.seqs))+1):end);


% testdata=struct;
% testdata.seqs=gooddatarand.seqs((round(0.9*length(gooddatarand.seqs))+1):end);
% testdata.VYBmus=gooddatarand.VYBmus((round(0.9*length(gooddatarand.seqs))+1):end);
% testdata.bincounts=gooddatarand.bincounts((round(0.9*length(gooddatarand.seqs))+1):end,:);
% 


data=[validdata  traindata];
outputfilename={'valid_enriched_rbzcontinuous_nomemorize_wstruct.mat','train_enriched_rbzcontinuous_nomemorize_wstruct.mat'};
%%
for k=1:length(data)
% plot enrichment versus mu
setfig('enrichment mu correlation');clf
plot(data(k).VYBmus,log2(sum(data(k).bincounts,2)),'.')
xlabel('\mu')
ylabel('coverage')

counts=log2(sum(data(k).bincounts,2));
muthresh=prctile(data(k).VYBmus,18);
mugood=data(k).VYBmus<muthresh;
enrichedcounts=counts;
enrichedcounts(mugood)=(enrichedcounts(mugood)+3); % arbitrarily enrich seqs according to their mu

hold on
plot(data(k).VYBmus,ceil(enrichedcounts),'.')
hold off


enrichedseqs={};
enrichedstructs={};

enrichedmus=[];

for i=1:length(data(k).seqs)
    ec=ceil(enrichedcounts(i));
    if data(k).VYBmus(i)<muthresh
    expandedseqs=repmat(data(k).seqs(i),ec,1);
    expandedstructs=repmat(data(k).structs(i),ec,1);
    expandedmus=repmat(data(k).VYBmus(i),ec,1);
    if length(enrichedseqs)>1
        enrichedseqs={enrichedseqs{:},expandedseqs{:}};
        enrichedstructs={enrichedstructs{:},expandedstructs{:}};
        enrichedmus=[enrichedmus(:);expandedmus];
    else
        enrichedseqs=expandedseqs;
        enrichedstructs=expandedstructs;
        enrichedmus=expandedmus;
    end
    else
        enrichedseqs{end+1}=data(k).seqs{i};
        enrichedstructs{end+1}=data(k).structs{i};
        enrichedmus(end+1)=data(k).VYBmus(i);
    end
end


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
structures=cell(length(enrichedstructs),1);
for i=1:length(enrichedseqs)
    s=enrichedseqs{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    seq{i}=s-'0';
    
    
    st=enrichedstructs{i};
    st=regexprep(st,'(','1');
    st=regexprep(st,'\.','2');
    st=regexprep(st,')','3');
    structures{i}=st-'0';
end

enriched.seqs=enrichedseqs;
enriched.mus=enrichedmus';
enriched.seqnum=seq;
enriched.structs=enrichedstructs;
enriched.structnum=structures;


totlen=113;
l1l2seqs=zeros(length(enriched.seqs),totlen);

for i=1:length(enriched.seqnum)
    l1l2seqs(i,1:length(enriched.seqnum{i}))=enriched.seqnum{i};
end

l1l2structs=zeros(length(enriched.structs),totlen);

for i=1:length(enriched.structnum)
    l1l2structs(i,1:length(enriched.structnum{i}))=enriched.structnum{i};
end
% convert to bit 
bitmat1=zeros(length(enriched.seqs), totlen,2,4);

for i=1:length(enriched.seqs)
    for l=1:totlen
        for j=1:4
            bitmat1(i,l,1,j)=l1l2seqs(i,l)==j;
        end
        for r=1:3
            bitmat1(i,l,2,r)=l1l2structs(i,l)==r;
        end
    end
end



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
%%
% v=struct;
% v.validxdata=bitmatrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))),:,:);
% v.validdata=switchesrand((round(0.8*length(bitmatrand(:,1,1)))+1):round(0.9*length(bitmatrand(:,1,1))));
% save('valid_enriched_rbzcontinuous_nomemorize.mat','v');
% %%
% tt=struct;
% tt.testxdata=bitmatrand((round(0.9*length(bitmatrand(:,1,1)))+1):end,:,:);
% tt.testdata=switchesrand((round(0.9*length(bitmatrand(:,1,1)))+1):end);
% save('test_enriched_rbzcontinuous.mat','tt');


%% Extract the loop sequences

loop1={};
loop2={};
loop1len=[];
loop2len=[];
for i=1:length(enriched.seqs)
    o=regexp(enriched.seqs{i},'(GCTGTCACCGGA)([A|C|T|G]*)(TCCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)','tokens');
    
    try
        loop1{end+1}=o{1}{2};
        loop2{end+1}=o{1}{4};
        loop1len(end+1)=length(o{1}{2});
        loop2len(end+1)=length(o{1}{4});
    catch
        loop1{end+1}='';
        loop2{end+1}='';
        loop1len(end+1)=0;
        loop2len(end+1)=0;
    end
end


enriched.loop1=loop1;
enriched.loop2=loop2;
enriched.loop1len=loop1len;
enriched.loop2len=loop2len;

%% go through all sequences and convert to matrix of 1 2 3 4
loop1seq=cell(length(enriched.seqs),1);
loop2seq=loop1seq;
for i=1:length(enriched.seqs)
    s=enriched.loop1{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop1seq{i}=s-'0';
    
    s=enriched.loop2{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop2seq{i}=s-'0';
    
end

enriched.loop1seqnum=loop1seq;
enriched.loop2seqnum=loop2seq;





%% pull out seqs with specified loop lengths
l1len=[5 6 7 8];
l2len=[15:25 30 40 50 60];
endslen=[3 3];

%% Analyses

muth=enriched.mus<prctile(enriched.mus,100);
[l1l2seqs,l1l2]=findSeqs(enriched,l1len,l2len,endslen,muth');

% what is the distribution of mus
findParameters(enriched.mus(l1l2),l1l2seqs);
%%
% find mutual sequence contribution - pairwise mu
[pairwisemu,pairwisesigma]=findPairwiseMu(enriched.mus(l1l2),l1l2seqs,[],1);

%%
% find entropy
ent=findEnt(l1l2seqs);

% find mutual information
[mI,pJoint,allI]=findMutualInformation(l1l2seqs);



%%









%% END




