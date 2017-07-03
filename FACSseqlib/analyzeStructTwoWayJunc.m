addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
g=load('~/Documents/YFS/YFSI_gooddataWstruct.mat');
g=g.gooddata;

%%
% twoway='(^\(*\.+)(\(+\.*\)+)(\.+\)*$)'; % includes variable stem length
% need a way to also include the nested loops, but exclude 3+ way juncs

twoway='(^\.+)(\(+.*\)+)(\.+$)';
twowaysidel='(^\.+)(\(+.*\)+$)';
twowaysider='(^\(+.*\)+)(\.+$)';

threeway='(^\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+$)';
threewaysidel='(^\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+$)';
threewaysider='(^\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+$)';

fourway='(^\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+$)';
fourwaysidel='(^\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+$)';
fourwaysider='(^\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+$)';

fiveway='(^\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+$)';
fivewaysidel='(^\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+$)';
fivewaysider='(^\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+)(\(+\.*\)+)(\.+$)';

loopstruct=g.loop2struct;

has2way=[find(~cellfun('isempty',regexp(loopstruct,twoway)) .* ...
               cellfun('isempty',regexp(loopstruct,threeway)) .* ...
               cellfun('isempty',regexp(loopstruct,fourway)) .* ...
               cellfun('isempty',regexp(loopstruct,fiveway))) ...
         find(~cellfun('isempty',regexp(loopstruct,twowaysidel)) .* ...
               cellfun('isempty',regexp(loopstruct,threewaysidel)) .* ...
               cellfun('isempty',regexp(loopstruct,fourwaysidel)) .* ...
               cellfun('isempty',regexp(loopstruct,fivewaysidel))) ...
         find(~cellfun('isempty',regexp(loopstruct,twowaysider)) .* ...
               cellfun('isempty',regexp(loopstruct,threewaysider)) .* ...
               cellfun('isempty',regexp(loopstruct,fourwaysider)) .* ...
               cellfun('isempty',regexp(loopstruct,fivewaysider)))];
struct2way=g.ensembstruct(find(has2way));

seq2way={};
loop1twoway={};
loop2twoway={};
loop1structtwoway={};
loop2structtwoway={};
for i=1:length(has2way)
    seq2way{end+1}=g.seqs{has2way(i)};
    loop1twoway{end+1}=g.loop1{has2way(i)};
    loop2twoway{end+1}=g.loop2{has2way(i)};
    loop1structtwoway{end+1}=g.loop1struct{has2way(i)};
    loop2structtwoway{end+1}=g.loop2struct{has2way(i)};
end
mus2way=g.mus(find(has2way));
length(seq2way)

setfig('twowayjunc mus');clf
hist(mus2way)

%% look for size of bulge
bulgeleftsize=[];
bulgerightsize=[];
for i=1:length(loop2structtwoway)
%     o=regexp(loop2structtwoway{i},'(^\(*)(\.+)(\(+.*\)+)(\.+)(\)*$)','tokens');
    o=regexp(loop2structtwoway{i},'(^\.*)(\(+\.*\)+)(\.*$)','tokens');
    try
        bulgeleftsize(end+1)=length(o{1}{1});
        
%         bulgeleftsize(end+1)=length(o{1}{2});
    catch
        bulgeleftsize(end+1)=0;
    end
    try
        bulgerightsize(end+1)=length(o{1}{3});
%         bulgerightsize(end+1)=length(o{1}{4});
    catch
        bulgerightsize(end+1)=0;
    end
end
%%
setfig('bulgesize');clf
subplot(1,3,1)
[n,c]=hist(bulgeleftsize,max(bulgeleftsize)-min(bulgeleftsize));
area(c+0.5,n)
set(gca,'yscale','log')
title('bulge left non-base paired')

subplot(1,3,2)
[n,c]=hist(bulgerightsize,max(bulgerightsize)-min(bulgerightsize));
area(c+0.5,n)
set(gca,'yscale','log')
title('bulge right non-base paired')

subplot(1,3,3)
[n,c]=hist(bulgeleftsize+bulgerightsize,max(bulgeleftsize+bulgerightsize)-min(bulgeleftsize+bulgerightsize));
area(c+0.5,n)
set(gca,'yscale','log')
title('bulge left+right non-base paired')

%% how many (a)symmetrical bulges are there?
buldgediff=bulgeleftsize-bulgerightsize;
setfig('symmetrical?');clf
[n,c]=hist(buldgediff,max(buldgediff)-min(buldgediff));
area(c+0.5,n)
set(gca,'yscale','log')
title('bulge left-right difference')
%%
symb=find(buldgediff==0);
symbl=bulgeleftsize(symb);
symbr=bulgerightsize(symb);

setfig('symmetrical bulgesize');clf
[n,c]=hist(symbl,max(symbl)-min(symbl));
area(c-0.5,n)
set(gca,'yscale','log')
title('symmetrical bulge lengths')

setfig('symmetrical bulgesize mus');clf

symbmus=mus2way(symb);
sl=[];
for i=min(symbl):max(symbl)
    sl(end+1)=mean(symbmus(find(symbl==i)));
end
hold on
plot(symbl,symbmus,'o','linewidth',2,'MarkerSize',5)
plot(min(symbl):max(symbl),sl,'-','linewidth',2)
hold off

legend('individual \mu values','mean \mu values','location','best')
xlabel('length of symmetrical bulge')
ylabel('\mu')


%% Note: Both theophylline and FA have symmetrical bulges, neo is off to 
% the side and tetracyline is a complicated 3 way junction
sym1bulge=find(symbr==1);
setfig('sym1bulge mu');clf
[n,c]=hist(mus2way(symb(sym1bulge)));
area(c-0.5,n)
title('symmetrical 1-base bulge \mu s')

loop2_1bulge={};
for i=1:length(sym1bulge)
    loop2_1bulge{end+1}=loop2structtwoway{symb(sym1bulge(i))};
end

sym2bulge=find(symbr==2);
setfig('sym2bulge mu');clf
[n,c]=hist(mus2way(symb(sym2bulge)));
area(c-0.5,n)
title('symmetrical 2-base bulge \mu s')

loop2_2bulge={};
for i=1:length(sym2bulge)
    loop2_2bulge{end+1}=loop2structtwoway{symb(sym2bulge(i))};
end

sym3bulge=find(symbr==3);
setfig('sym3bulge mu');clf
[n,c]=hist(mus2way(symb(sym3bulge)));
area(c-0.5,n)
title('symmetrical 3-base bulge \mu s')

loop2_3bulge={};
for i=1:length(sym3bulge)
    loop2_3bulge{end+1}=loop2structtwoway{symb(sym3bulge(i))};
end


sym5bulge=find(symbr==5);
setfig('sym5bulge mu');clf
[n,c]=hist(mus2way(symb(sym5bulge)));
area(c-0.5,n)
title('symmetrical 5-base bulge \mu s')

loop2_5bulge={};
for i=1:length(sym5bulge)
    loop2_5bulge{end+1}=loop2structtwoway{symb(sym5bulge(i))};
end

%% END



