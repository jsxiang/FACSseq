function [outstructs, ind, totallooplength]=findStructs(gooddata,loop1,loop2,endslen,percentile)
% objective is to give as input loop sizes, and output the corresponding
% sequences and indices for gooddata.m

if length(loop1) >1 && length(loop2) >1
    totallooplength=sum(endslen)*2;
    l1ind=zeros(1,length(gooddata.seqs));
    for i=1:length(loop1)
        l1ind=l1ind| gooddata.loop1len==loop1(i);
    end
    l2ind=zeros(1,length(gooddata.seqs));
    for i=1:length(loop2)
        l2ind=l2ind| gooddata.loop2len==loop2(i);
    end
    ind=find(l1ind.*l2ind.*percentile);
    outstructs=zeros(length(ind),totallooplength);

    for i=1:length(ind)
        l1=gooddata.loop1structnum{ind(i)};
        l1start=l1(1:endslen(1));
        l1end=l1((end-endslen(2)+1):end);
        
        l2=gooddata.loop2structnum{ind(i)};
        l2start=l2(1:endslen(1));
        l2end=l2((end-endslen(2)+1):end);
        
%         outseqs(i,1:endslen(1))=l1start;
%         outseqs(i,(endslen(1)+1):(endslen(1)+endslen(2)))=l1end;
% 
%         outseqs(i,(loop1+1):(loop1+endslen(1)))=l2start;
%         outseqs(i,(loop1+endslen(1)+1):(loop1+endslen(1)+endslen(2)))=l2end;
        outstructs(i,:)=[l1start l1end l2start l2end];
    end
elseif length(loop2)>1 || loop2(1)>10
    totallooplength=loop1+sum(endslen);
    l1ind=zeros(1,length(gooddata.seqs));
    for i=1:length(loop1)
        l1ind=l1ind| gooddata.loop1len==loop1(i);
    end
    
    l2ind=zeros(1,length(gooddata.seqs));
    for i=1:length(loop2)
        l2ind=l2ind| gooddata.loop2len==loop2(i);
    end
    
    ind=find(l1ind.*l2ind.*percentile);
    outstructs=zeros(length(ind),totallooplength);
    size(ind)
    for i=1:length(ind)
        l2=gooddata.loop2structnum{ind(i)};
        l2start=l2(1:endslen(1));
        l2end=l2((end-endslen(2)+1):end);

        outstructs(i,(loop1+1):(loop1+endslen(1)))=l2start;
        outstructs(i,(loop1+endslen(1)+1):(loop1+endslen(1)+endslen(2)))=l2end;
%         fprintf('index = %0.0f, loop1len = %0.0f\n',ind(i),loop1);

        outstructs(i,1:loop1)=gooddata.loop1structnum{ind(i)};
    end
elseif length(loop1)>1 || loop1(1)>10

    totallooplength=sum(endslen)+loop2;
    l1ind=zeros(1,length(gooddata.seqs));   
    for i=1:length(loop1)
        l1ind=l1ind| gooddata.loop1len==loop1(i);
    end

    l2ind=zeros(1,length(gooddata.seqs));
    for i=1:length(loop2)
        l2ind=l2ind| gooddata.loop2len==loop2(i);
    end

    ind=find(l1ind.*l2ind.*percentile);
    outstructs=zeros(length(ind),totallooplength);
    for i=1:length(ind)
        l1=gooddata.loop1structnum{ind(i)};
        l1start=l1(1:endslen(1));
        l1end=l1((end-endslen(2)+1):end);

        outstructs(i,1:endslen(1))=l1start;
        outstructs(i,(endslen(1)+1):(endslen(1)+endslen(2)))=l1end;

        outstructs(i,(sum(endslen)+1):(sum(endslen)+loop2))=gooddata.loop2structnum{ind(i)};
    end

end