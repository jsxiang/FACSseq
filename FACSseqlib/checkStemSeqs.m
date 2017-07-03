addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS

%%
g=load('YFSI_gooddata.mat');
gooddata=g.gooddata;

%% when the sequence contains the beginning and end of theoaptamer 

theoAAG='ATACCAGCATCGTCTTGATGCCCTTGGAAG';
theoCAG='ATACCAGCATCGTCTTGATGCCCTTGGCAG';
neo='GCTTGTCCTTTAATGGTCC';
tet='AAAACATACCAGATTTCGATCTGGAGAGGTGAAGAATTCGACCACCT';
FA='GTGCTTGGTACGTTATATTCAGC';
aptamer=theoCAG;
aptamername='theoCAG';

%%
matchstring1='^G.*GC$';
matchstring2='^G.*GC$';


aptONl1=zeros(1,length(gooddata.seqs));
aptONl2=zeros(1,length(gooddata.seqs));
for i=1:length(gooddata.seqs)
    if gooddata.loop1len(i)>14 & regexp(gooddata.loop1{i},matchstring1)
        aptONl1(i)=1;
    elseif gooddata.loop2len(i)>14 & regexp(gooddata.loop2{i},matchstring2)
        aptONl2(i)=1;
    end
end
aptONl1=find(aptONl1>0);
aptONl2=find(aptONl2>0);

setfig('apt on Loop1');clf

aptstemONl1.mus=[];
aptstemONl1.seqs={};
aptstemONl1.l1={};
aptstemONl1.l2={};
for i=1:length(aptONl1)
aptstemONl1.seqs{end+1}=gooddata.seqs{aptONl1(i)};
aptstemONl1.mus(end+1)=gooddata.mus(aptONl1(i));
aptstemONl1.l1{end+1}=gooddata.loop1{aptONl1(i)};
aptstemONl1.l2{end+1}=gooddata.loop2{aptONl1(i)};
subplot(length(aptONl1),2,2*i-1)
bar(gooddata.bincounts(aptONl1(i),1:12))
subplot(length(aptONl1),2,2*i)
bar(gooddata.bincounts(aptONl1(i),13:24))

end

setfig('apt on Loop2');clf

aptstemONl2.mus=[];
aptstemONl2.seqs={};
aptstemONl2.l1={};
aptstemONl2.l2={};
for i=1:length(aptONl2)
aptstemONl2.seqs{end+1}=gooddata.seqs{aptONl2(i)};
aptstemONl2.mus(end+1)=gooddata.mus(aptONl2(i));
aptstemONl2.l1{end+1}=gooddata.loop1{aptONl2(i)};
aptstemONl2.l2{end+1}=gooddata.loop2{aptONl2(i)};
subplot(length(aptONl2),2,2*i-1)
bar(gooddata.bincounts(aptONl2(i),1:12))
subplot(length(aptONl2),2,2*i)
bar(gooddata.bincounts(aptONl2(i),13:24))

end
