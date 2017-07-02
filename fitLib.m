clear
addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS

gcstruct=load('MFSI_goodcoveragerep2.mat');
gc=gcstruct.goodcoverage;

FA_apt='GC[A|C|G|T]TGGTACGTTATATTC[A|G]G';
theo_apt='ATACCAGCATCGTCTTGATGCCCTTGGAAG';

theo=struct; % initiate empty struct
FA=struct;
for i=1:length(gc)
        matchresult=regexp(gc(i).seq,theo_apt, 'once');
        if ~isempty(matchresult{1})
            theo(end+1).seq=gc(i).seq;
            theo(end).cond=gc(i).cond;
        end
        
        matchresult=regexp(gc(i).seq,FA_apt, 'once');
        if ~isempty(matchresult{1})
            FA(end+1).seq=gc(i).seq;
            FA(end).cond=gc(i).cond;
        end
        
end
%%
fid=fopen('../libreads/MFSI_libnoBC_MoreThanOnce.txt');
seqs=textscan(fid,'%s');
fclose(fid);

dat=load('MFSI.mat'); % this is loaded in as a struct

x=load('barcodeglobalcount.txt');
scalex=(x(1:40,1)./x(1:40,3))';
scalexmat=repmat(scalex,length(seqs{:}),1);
scaleddata=dat.MFSIdat(:,1:40).*scalexmat;
dat.seqs=seqs{:};

dat.cond(1).scaled=scaleddata(:,1:10);
dat.cond(2).scaled=scaleddata(:,11:20);
dat.cond(3).scaled=scaleddata(:,21:30);
dat.cond(4).scaled=scaleddata(:,31:40);

save('MFSI_scaled.mat','dat')

%%
theoAll=struct;
theo_apt='ATACCAGCATCGTCTTGATGCCCTTGGAAG';

for i=2:length(dat.seqs)
    matchresult=regexp(dat.seqs(i),theo_apt, 'once');
        if ~isempty(matchresult{1})
            theoAll(end+1).seq=dat.seqs(i);
            for k=1:4
                theoAll(end).cond(k).data=dat.cond(k).scaled(i,:);
            end
        end
end

theomus_All=zeros(length(theoAll),4);
theosigmas=theomus_All;
for j=1:length(theoAll)
    for k=1:length(theoAll(j).cond)
        scaledbincounts=theoAll(j).cond(k).data;
        edges=1:10;
        [m,s]=normfit(edges,[],[],scaledbincounts);
        theoAll(j).cond(k).mu=m;
        theoAll(j).cond(k).sigma=s;
        theomus_All(j,k)=m;
        theosigmas(j,k)=s;         

    end
end
theomus_All(isnan(theomus_All))=0;
save('theoAll.mat','theoAll')
%% grab data and fit directly 
% going through all this data only took 3.123 seconds

FAmus=zeros(length(FA),4);
FAsigmas=FAmus;
for j=1:length(FA)
    for k=1:length(FA(j).cond)
        scaledbincounts=FA(j).cond(k).data;
        edges=1:10;
        
%         if sum(scaledbincounts(1)) ==0
%             continue
%         else
        [m,s]=normfit(edges,[],[],scaledbincounts);
        FA(j).cond(k).mu=m;
        FA(j).cond(k).sigma=s;
        FAmus(j,k)=m;
        FAsigmas(j,k)=s;
%         end

    end
end

sum(isnan(FAmus));
FAmus(isnan(FAmus))=0;

save('FA2.mat','FA')
%% grab data and fit directly 
% going through all this data only took 3.123 seconds

theomus=zeros(length(theo),4);
theosigmas=theomus;
for j=1:length(theo)
    for k=1:length(theo(j).cond)
        scaledbincounts=theo(j).cond(k).data;
        edges=1:10;
        
        
        [m,s]=normfit(edges,[],[],scaledbincounts);
        theo(j).cond(k).mu=m;
        theo(j).cond(k).sigma=s;
        theomus(j,k)=m;
        theosigmas(j,k)=s;         

    end
end
theomus(isnan(theomus))=0;
save('theo2.mat','theo')

%% comparing replicates
setfig('compare replicates FA');clf
x=linspace(min(edges),max(edges));
y=x;

allmus=FAmus;
subplot(2,2,1)
dscatter(allmus(:,1),allmus(:,3),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('2-')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)


subplot(2,2,2)
dscatter(allmus(:,2),allmus(:,4),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('1+')
ylabel('2+')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)


subplot(2,2,3)
dscatter(allmus(:,1),allmus(:,2),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('1+')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

subplot(2,2,4)
dscatter(allmus(:,3),allmus(:,4),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('2-')
ylabel('2+')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

%% comparing replicates
setfig('compare replicates theo');clf
x=linspace(min(edges),max(edges));
y=x;

allmus=theomus;
subplot(2,2,1)
dscatter(allmus(:,1),allmus(:,3),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('2-')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)


subplot(2,2,2)
dscatter(allmus(:,2),allmus(:,4),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('1+')
ylabel('2+')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)


subplot(2,2,3)
dscatter(allmus(:,1),allmus(:,2),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('1-')
ylabel('1+')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

subplot(2,2,4)
dscatter(allmus(:,3),allmus(:,4),'marker','o','BINS',[500 500])
hold on
plot(x,y,'k:','linewidth',1.5)
hold off
xlabel('2-')
ylabel('2+')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

