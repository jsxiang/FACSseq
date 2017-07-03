addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS

%%
g=load('YFSI_gooddata.mat');
gooddata=g.gooddata;
%%
% n8=load('/Volumes/smolke-lab$/Joy/NGS data backup/NGS8/ngs8.mat');
n10=load('/Volumes/smolke-lab$/Joy/NGS data backup/NGS10/ngs10vivo.mat');
%%
% n8ratios=n8.invivo.ratio;
% n8seqs={};
% for i=1:length(n8.tmpl.elib)
%     n8seqs{end+1}=n8.tmpl.elib(i).seq;
% end

%% This takes a long time to load

n10minus=n10.invivo(5).ratio;
n10plus=n10.invivo(6).ratio;
n10seqs={};
n10l1={};
n10l1len=[];

n10l2={};
n10l2len=[];
for i=1:length(n10.tmpl.elib)
    n10seqs{end+1}=n10.tmpl.elib(i).seq;
    n10l1{end+1}=n10.tmpl.elib(i).s1;
    n10l1len(end+1)=length(n10.tmpl.elib(i).s1);
    n10l2{end+1}=n10.tmpl.elib(i).s2;
    n10l2len(end+1)=length(n10.tmpl.elib(i).s2);
    
end


%% n10 
setfig('lengths dist n10');clf
m=zeros(85,85);
mc=m;
zerolengthloops=[];
r=round(rand(1,100000)*2030967);

for k=1:length(r)
    if n10l1len(r(k))==0||n10l2len(r(k))==0
        zerolengthloops(end+1)=k;
    else
        m(n10l1len(r(k)),n10l1len(r(k)),end+1)=n10minus(k);
        mc(n10l1len(r(k)),n10l2len(r(k)))=mc(n10l1len(r(k)),n10l2len(r(k)))+1;
    end
end
imagesc(mc)
 


%% when the sequence contains the beginning and end of theoaptamer 

theoAAG='ATACCAGCATCGTCTTGATGCCCTTGGAAG';
theoCAG='ATACCAGCATCGTCTTGATGCCCTTGGCAG';
neo='GCTTGTCCTTTAATGGTCC';
tet='AAAACATACCAGATTTCGATCTGGAGAGGTGAAGAATTCGACCACCT';
aptamer=theoAAG;
aptamername='theoAAG';

n10l1long = cellfun(@(s) ~isempty(strfind(aptamer, s)), n10l1);
n10l2long = cellfun(@(s) ~isempty(strfind(aptamer, s)), n10l2);

n10l1longl2seq=n10l2(find(n10l1long));
n10l1longmu=containers.Map(n10l1longl2seq,n10minus(find(n10l1long)));
n10l2longl1seq=n10l1(find(n10l2long));
n10l2longmu=containers.Map(n10l2longl1seq,n10minus(find(n10l2long)));
%%
matchstring1='^A.*AAG$';
matchstring2='^A.*GAAG$';


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

%%
matchedn10=[];
matchedYFSI=[];
for i=1:length(aptstemONl1.l2)
    try
        matchedn10(end+1)=n10l1longmu(aptstemONl1.l2{i});
        matchedYFSI(end+1)=aptstemONl1.mus(i);
        
    catch
    end
end

for i=1:length(aptstemONl2.l1)
    try
        matchedn10(end+1)=n10l1longmu(aptstemONl2.l1{i});
        matchedYFSI(end+1)=aptstemONl2.mus(i);
    catch
    end
end
setfig('compare n10 vs YFSI');clf
try
    plot(matchedn10,matchedYFSI,'o')
    mdl = fitlm(matchedn10,matchedYFSI);
    rsq=mdl.Rsquared.Adjusted;
    t=sprintf('R^2 = %0.2f',rsq);
    text(0.02,10,t)
    xlabel('NGS10')
    ylabel('YFSI')
    hold on
    x=linspace(min(matchedn10),max(matchedn10));
    y=x;
    plot(x,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*x,'k:','linewidth',1.5)
    hold off
    tt=sprintf('%s matched with loop1=%s and loop2=%s',aptamername,matchstring1,matchstring2);
    title(tt,'interpreter','none')
catch
    fprintf('No valid sequences / mu found! \n')
end

%% END

