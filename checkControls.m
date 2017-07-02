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

%% graded ribozyme controls
ctrls={'strsv','TGTGCTT','GTGA';
'g811','TCAG','GTGA';
'g814','TGTT','ATAA';
'g816','TGCT','GTGA';
'g824','TGTT','ACTA';
'g830','GGAG','AAAT';
'g832','GGCT','AGCT';
'g833','TACT','CAGA';
'g836','CTTT','CAGA';
'g862','TGCA','CGCG'};

%%
CI=[];
for i=1:length(ctrls)
    matchstring1=ctrls{i,2};
    matchstring2=ctrls{i,3};

    ci1 = cellfun(@(s) ~isempty(strfind(matchstring1, s)), gooddata.loop1);
    ci2 = cellfun(@(s) ~isempty(strfind(matchstring2, s)), gooddata.loop2);

    ci=ci1&ci2;
    try
        CI(end+1)=find(ci>0);
    catch
        if isempty(find(ci>0))
            CI(end+1)=nan;
        end
    end
end
%%
valided={'YI14379','YI41145','YI40685','YI21508','YI9874','YI7971','YI4263','YI15657','YI39423','YI38802','YI26640','YI3102'};
gmus=containers.Map(gooddata.index,gooddata.mus);
%%
gseqs=containers.Map(gooddata.index,1:length(gooddata.seqs));
%%
gocounts=containers.Map(gooddata.index,1:length(gooddata.origcounts));
ycounts=zeros(length(valided),length(gooddata.origcounts(1,:)));
yseqs={};
ymus=[]
for i=1:length(valided)
    try
        ycounts(i,:)=gooddata.origcounts(gocounts(valided{i}),:);
    end
    yseqs{end+1}=gooddata.seqs{gseqs(valided{i})};
    ymus(end+1)=gmus(valided{i});
end

totalcounts=sum(ycounts,2)
for i=1:length(totalcounts)
    fprintf('%0.0f %s\n',totalcounts(i),yseqs{i})
end

%% from grepbincounts.xlsx

bc=[320	328	105	116	67	51	19	6	9	4	0	1	;
0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	30	45	0	0	1	0	2	0	0	0	;
0	2	25	101	3	0	5	0	0	0	0	0	;
0	0	11	39	65	50	32	31	0	0	0	0	;
0	0	0	6	76	77	69	55	25	1	0	0	;
0	2	6	10	29	82	29	55	249	150	49	3	;
1	11	6	7	3	10	7	9	46	157	135	14	;
0	2	1	9	0	10	2	10	11	46	192	103	;
0	2	0	0	0	2	2	0	3	0	10	50	;
0	0	1	0	0	0	3	0	0	0	0	1	;
169	140	47	40	16	15	5	2	1	2	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	7	0	0	0	0	0	0	0	0	0	;
0	0	26	17	0	0	0	0	0	0	0	0	;
0	2	51	99	0	2	0	0	0	0	0	0	;
2	0	28	94	62	52	88	33	0	0	0	0	;
0	0	0	5	74	97	87	131	91	4	2	0	;
0	0	3	2	9	23	13	32	71	40	5	1	;
0	7	2	1	0	4	6	1	21	73	69	4	;
0	3	0	0	0	6	1	9	0	23	70	57	;
0	0	1	0	1	0	0	0	0	3	17	89	;
1	2	0	4	0	0	1	0	0	7	4	3	]';

ycounts==bc;


% setfig('apt on Loop1');clf
% 
% aptstemONl1.mus=[];
% aptstemONl1.seqs={};
% aptstemONl1.l1={};
% aptstemONl1.l2={};
% for i=1:length(aptONl1)
% aptstemONl1.seqs{end+1}=gooddata.seqs{aptONl1(i)};
% aptstemONl1.mus(end+1)=gooddata.mus(aptONl1(i));
% aptstemONl1.l1{end+1}=gooddata.loop1{aptONl1(i)};
% aptstemONl1.l2{end+1}=gooddata.loop2{aptONl1(i)};
% subplot(length(aptONl1),2,2*i-1)
% bar(gooddata.bincounts(aptONl1(i),1:12))
% subplot(length(aptONl1),2,2*i)
% bar(gooddata.bincounts(aptONl1(i),13:24))
% 
% end
% 
% setfig('apt on Loop2');clf
% 
% aptstemONl2.mus=[];
% aptstemONl2.seqs={};
% aptstemONl2.l1={};
% aptstemONl2.l2={};
% for i=1:length(aptONl2)
% aptstemONl2.seqs{end+1}=gooddata.seqs{aptONl2(i)};
% aptstemONl2.mus(end+1)=gooddata.mus(aptONl2(i));
% aptstemONl2.l1{end+1}=gooddata.loop1{aptONl2(i)};
% aptstemONl2.l2{end+1}=gooddata.loop2{aptONl2(i)};
% subplot(length(aptONl2),2,2*i-1)
% bar(gooddata.bincounts(aptONl2(i),1:12))
% subplot(length(aptONl2),2,2*i)
% bar(gooddata.bincounts(aptONl2(i),13:24))
% 
% end

