function gooddata=findSeqnum(gooddata)

for j=1:length(gooddata)
    loop1seq=cell(length(gooddata(j).seqs),1);
    loop2seq=loop1seq;
for i=1:length(gooddata(j).seqs)
    s=gooddata(j).loop1{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop1seq{i}=s-'0';
    
    s=gooddata(j).loop2{i};
    s=regexprep(s,'A','1');
    s=regexprep(s,'T','2');
    s=regexprep(s,'C','3');
    s=regexprep(s,'G','4');
    % this is a super cool way of converting a string of numbers to a matrix
    loop2seq{i}=s-'0';
end


gooddata(j).loop1seqnum=loop1seq;
gooddata(j).loop2seqnum=loop2seq;
end
end