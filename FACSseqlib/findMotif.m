function data=findMotif(inputdata,motifs,motifname)
    data=struct;
    for i=1:length(motifs)
    s=regexp(inputdata.seqs,motifs(i),'tokens');
    hasmotif=~cellfun('isempty',s);
    data(i).mus=inputdata.mus(hasmotif,:);
    
    data(i).sigma=inputdata.mus(hasmotif,:);
    data(i).rep1ind=inputdata.rep1ind;
    data(i).rep2ind=inputdata.rep2ind;
    try
    data(i).rep3ind=inputdata.rep3ind;
    data(i).rep4ind=inputdata.rep4ind;
    end
        data(i).origbincounts=inputdata.origbincounts(hasmotif,:);

    try
        data(i).VYBmus=inputdata.VYBmus(hasmotif,:);
        data(i).ttestdiff=inputdata.ttestdiff(hasmotif,:);
    catch
    end
    data(i).seqs=inputdata.seqs(hasmotif);
    data(i).motifname=motifname(i);
    data(i).hasmotif=hasmotif;
    data(i).loop1={};
    data(i).loop2={};
    data(i).loop1len=[];
    data(i).loop2len=[];
    re=s(hasmotif);

        for k=1:length(data(i).seqs)
            try
            data(i).loop1{end+1}=re{k}{1}{2};
            data(i).loop2{end+1}=re{k}{1}{4};
            data(i).loop1len(end+1)=length(data(i).loop1{k});
            data(i).loop2len(end+1)=length(data(i).loop2{k});
            catch
            data(i).loop1{end+1}='';
            data(i).loop2{end+1}='';
            data(i).loop1len(end+1)=nan;
            data(i).loop2len(end+1)=nan;
            end

        end
    end
end