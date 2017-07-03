clear
addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/FACS/
% addpath ~/Documents/MATLAB/BREWER/
%% load data that are deemed "good", ones that have enough coverage
g=load('YFSI_gooddata.mat');
gooddata=g.gooddata;


%% go through all sequences and grab MFE from RNAstructure
ensembct={};
ensembstruct={};
ensembstructfrac=[];
loop1struct={};
loop2struct={};
fraccleavers=[];
for i=1:length(gooddata.seqs)
    seed=randi(10000);
    cnt=100;
    seq=gooddata.seqs{i};
    
    fid=fopen('_seq.fasta','w');
    fprintf(fid,'>seq\n%s',seq);
    fclose(fid);
    [s,r]=system(sprintf('export DATAPATH=~/Documents/bioinformatics/RNAStructure/data_tables; ~/Documents/bioinformatics/RNAStructure/exe/stochastic --sequence -s %d -e %d %s %s', seed, cnt, '_seq.fasta', '_seqfolded.ct'));
    % convert dot ct file to dot bracket
    ctready=1;
    while ctready
        try
        seqct=ctread('_seqfolded.ct');
        ctready=0;
        catch
        pause(0.1);
        end
    end
    ensembct{end+1}={seqct.folding};
    
    [c,ia,ic]=unique({seqct.folding});
    % the following four lines is to get the occurrences of each unique
    % sequence
    y=sort(ic);
    p = find([numel(y);diff(y);numel(y)]);
    values = y(p(1:end-1));
    instances = diff(p);
    ensembstruct{end+1}=c{instances==max(instances)};
    ensembstructfrac(end+1)=max(instances)/cnt;
    % remember to grab the structure of the loop sequences
    try
        o=regexp(ensembstruct{i},'(\(\(\(\(\(\.\(\(\(\(\(\()(.*)(\)\)\)\)\)\)\.\.\.\.\.\.\.\(\(\(\()(.*)(\)\)\)\)\.\.\.\)\)\)\)\))','tokens');
        loop1struct{end+1}=o{1}{2};
        loop2struct{end+1}=o{1}{4};
    catch
        loop1struct{end+1}='';
        loop2struct{end+1}='';
    end

    % define the active ribozyme scaffold
    rbz1='^(((((\.((((((';
    rbz2='))))\.\.\.)))))$';
    seqRE1=regexp({seqct.folding},rbz1);
    seqRE2=regexp({seqct.folding},rbz2);
    fraccleavers(end+1)=sum(~cellfun('isempty',seqRE1) .* ~cellfun('isempty',seqRE2))/cnt;
    system(sprintf('rm _seq*'));
    
    if rem(i,1000)==999

        gooddata.ensembct=ensembct;
        gooddata.ensembstruct=ensembstruct;
        gooddata.ensembstructfrac=ensembstructfrac;
        gooddata.loop1struct=loop1struct;
        gooddata.loop2struct=loop2struct;
        gooddata.fraccleavers=fraccleavers;

        save('YFSI_gooddata.mat','gooddata');
    end
    
end


%%

gooddata.ensembct=ensembct;
gooddata.ensembstruct=ensembstruct;
gooddata.ensembstructfrac=ensembstructfrac;
gooddata.loop1struct=loop1struct;
gooddata.loop2struct=loop2struct;
gooddata.fraccleavers=fraccleavers;
save('YFSI_gooddata.mat','gooddata');



%% END
