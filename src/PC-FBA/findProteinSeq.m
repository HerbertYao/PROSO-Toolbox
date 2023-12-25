function seq = findProteinSeq(geneID, fasta, defaultLen)

% The function returns a protein sequence from fasta file. If the ID is not
% found, return an average-lengthed 'ZZZZZZZZ' seq
% 
% USAGE:
%   
%   seq = findProteinSeq(geneID, fasta, defaultLen)
% 
% INPUTS:
%   geneID:     Query gene id to find in the fasta header
%   fasta:      Protein fasta file in the form of N*1 struct
%   defaultLen: The default aa length if the query is not found
% 
% OUTPUTS:
%   seq: Returned query sequence
% 
% EXAMPLE:
% 
%   seq = findProteinSeq('b0001', ecoli_fasta, 200)
%  
% .. AUTHOR: - Herbert Yao, Dec 2023
% 

found = false;

% First try complete match
for i = 1:length(fasta)
    if strcmp(fasta(i).Header,geneID)
        seq = fasta(i).Sequence;
        found = true;
        break;
    end
end

% Then try contains
if ~found
    for i = 1:length(fasta)
        if contains(fasta(i).Header,geneID)
            seq = fasta(i).Sequence;
            found = true;
            break;
        end
    end
end

% If still not found, calculate the average length and return
if ~found
    seq(1:round(defaultLen)) = 'Z';
    fprintf('Not found gene: %s\n',geneID);
end

end
    