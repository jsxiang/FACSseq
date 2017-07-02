addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS
addpath ~/Documents/MATLAB/FACSseq/
clear 
c=load('commondata.mat');
commondata=c.commondata;

%%
switchref=1; %1 is theo, 2 is FA
%%
largeL1={'(GCTGTCACCGG)([A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G]*)(CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)'};

d=~cellfun('isempty',regexp(commondata(switchref).seqs,largeL1));

nonswitchingTheo=commondata(switchref).seqs(find((d & ~commondata(switchref).switches)));
theoRS(1).seqs=nonswitchingTheo;
theoRS(1).sigma=commondata(switchref).semR1(find((d & ~commondata(switchref).switches)));
theoRS(1).mus=commondata(switchref).yVYBmus1(find((d & ~commondata(switchref).switches)))';
theo(1)=findMotif(theoRS(1),largeL1,{'nonswitchesLargeL1'});


theosn(1)=findSeqnum(theo(1));
%% Analyses
i=1;
data=theosn;
l2len=[5];
l1len=[32]
endslen=[10 9];

muthlo=data(1).mus<prctile(max(data(2).mus),100); % filter
muthhi=data(2).mus<prctile(max(data(2).mus),100); % filter
% muth=~cellfun('isempty',regexp(gooddata.loop1,'^T[A|T|C|G]+'));
[l1l2seqslo,l1l2lo]=findSeqs(data(1),l1len,l2len,endslen,muthlo');
[l1l2seqshi,l1l2hi]=findSeqs(data(2),l1len,l2len,endslen,muthhi');
% seqlogo(data(i).seqs(l1l2lo))
% % what is the distribution of mus
% findParameters(gooddata.VYBmus(l1l2),l1l2seqs);
% find mutual sequence contribution - pairwise mu

set(gca,'fontsize',14)

% find entropy
[entlo,plo]=findEnt(l1l2seqslo,data(i).motifname);
[enthi,phi]=findEnt(l1l2seqshi,data(i).motifname);

[dent,dp]=findDeltaEnt(entlo,enthi,plo,phi,data(i).motifname);

%% find mutual information
[mIlo,pJointlo,allIlo]=findMutualInformation(l1l2seqslo,data(i).motifname);
%%
[mIhi,pJointhi,allIhi]=findMutualInformation(l1l2seqshi,data(i).motifname);
%%
[deltaMI,deltapJoint]=findDeltaMI(mIlo,mIhi,pJointlo,pJointhi,data(i).motifname);


%%

largeL1={'(GCTGTCACCGG)([A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G]*)(CCGGTCTGATGAGTCC)([A|C|T|G]*)(GGACGAAACAGC)'};

d=~cellfun('isempty',regexp(commondata(switchref).seqs,largeL1));




allcommonTheoswitches=commondata(switchref).seqs(find((d & commondata(switchref).switches)));
theoRS(2).seqs=allcommonTheoswitches;
theoRS(2).sigma=commondata(switchref).semR2(find((d & commondata(switchref).switches)));
theoRS(2).mus=commondata(switchref).yVYBmus1(find((d & commondata(switchref).switches)))';
theo(2)=findMotif(theoRS(2),largeL1,{'switchesLargeL1'});

theosn(2)=findSeqnum(theo(2));
%% Analyses
i=2;
data(i)=theosn(i);



l2len=[5];
l1len=[19:24]
endslen=[9 9];

muthlo=data(i).mus<prctile(data(i).mus,40); % filter
muthhi=data(i).mus>prctile(data(i).mus,60); % filter
% muth=~cellfun('isempty',regexp(gooddata.loop1,'^T[A|T|C|G]+'));
[l1l2seqslo,l1l2lo]=findSeqs(data(i),l1len,l2len,endslen,muthlo');
[l1l2seqshi,l1l2hi]=findSeqs(data(i),l1len,l2len,endslen,muthhi');
% seqlogo(data(i).seqs(l1l2lo))
% % what is the distribution of mus
% findParameters(gooddata.VYBmus(l1l2),l1l2seqs);
% find mutual sequence contribution - pairwise mu
[pairwisemulo,pairwisesigmalo]=findPairwiseMu(data(i).mus(l1l2lo),l1l2seqslo,10,1,data(i).motifname);
% pairwise50=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,50,[]);
% pairwise25=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,25,[]);
% pairwise10=findPairwiseMu(gooddata.VYBmus(l1l2),l1l2seqs,10,[]);

set(gca,'fontsize',14)

% find entropy

[entlo,plo]=findEnt(l1l2seqslo,data(i).motifname);
[enthi,phi]=findEnt(l1l2seqshi,data(i).motifname);

[dent,dp]=findDeltaEnt(entlo,enthi,plo,phi,data(i).motifname);


%% find mutual information
[mIlo,pJointlo,allIlo]=findMutualInformation(l1l2seqslo,data(i).motifname);
[mIhi,pJointhi,allIhi]=findMutualInformation(l1l2seqshi,data(i).motifname);
%%
[deltaMI,deltapJoint]=findDeltaMI(mIlo,mIhi,pJointlo,pJointhi,data(i).motifname);

%%
setfig('distribution of mus');clf
hold on
h1=histogram(theosn(1).mus);
h2=histogram(theosn(2).mus);
h1.Normalization='probability';
h2.Normalization='probability';
h1.NumBins=100;
h2.NumBins=40;


%%
f=fopen('nonswitchingTheo.fasta','w');
fprintf(f,'>seq\n%s\n',nonswitchingTheo{:});
fclose(f);

%%



f=fopen('allswitchingTheo.fasta','w');
for i=1:length(allcommonTheoswitches)
    fprintf(f,'>%s_%1.2f\n%s\n',theosn(2).loop2{i},theosn(2).mus(i),theosn(2).seqs{i});
end
fclose(f);












