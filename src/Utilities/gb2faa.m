function faa = gb2faa(gbFileName,ifDropNoSeqGene,faaFileName)

% Parse a protein FASTA file from GenBank file
% 
% USAGE:
% 
%   faa = gb2faa(gbFileName)
%   faa = gb2faa(gbFileName,ifDropNoSeqGene)
%   faa = gb2faa(gbFileName,ifDropNoSeqGene,faaFileName)
%   
% 
% INPUTS:
%   gbFileName: Genbank file name. May also include its path
% 
% OPTIONAL INPUTS: 
%   ifDropNoSeqGene: if genes without sequence recorded are passed. This is
%                    recommended because faa with empty seq may trigger
%                    errors downstream. Default: true
%   faaFileName:     Output faa file name. Faa file will not be created if
%                    this is left blank
% 
% OUTPUTS: 
%   faa: faa output format in N*1 cell. Odd entries are gene headings
%        (attributes) starting with '>', even entries are gene sequences. 
% 
% EXAMPLE:
% 
%   faa = gb2faa('GCF_000005845.2_ASM584v2_genomic.gbff.tar',true,'GCF_000005845.2_ASM584v2');
% 
% .. AUTHOR: - Herbert Yao, Dec 2023
% 

if ~exist('dropNoSeqGene','var')
    ifDropNoSeqGene = true;
end

gb = genbankread(gbFileName);

faa = {};

for i = 1:length(gb) % chromosome
    for j = 1:length(gb(i).CDS) % gene
        
%       record seq

        if ~isempty(gb(i).CDS(j).translation)
            seq = gb(i).CDS(j).translation;
        elseif ~ifDropNoSeqGene
            seq = '';
        else
            continue; % we suggest to pass entries with no seq
        end

%       record attributes from CDS

        attr = {};
        attrVal = {};

        txt = gb(i).CDS(j).text;
        [hei,~] = size(txt);

        for k = 1:hei % parse from gb

            if startsWith(txt(k,:),'/translation')
                break; % break when reach translation rows

            elseif k == 1 % location
                attr{end+1} = 'location';
                splt = erase(txt(k,:),'CDS');
                attrVal{end+1} = erase(splt,' ');

            elseif startsWith(txt(k,:),'/') % add new attribute
                try
                    if contains(txt(k,:),'"') % text-based
                        splt = split(txt(k,:),'"');
                        attrVal{end+1} = splt{2};
                        attr{end+1} = splt{1}(2:end-1);
                    else % number-based
                        splt = split(txt(k,:),'=');
                        attrVal{end+1} = erase(splt{2},' ');
                        attr{end+1} = splt{1}(2:end);
                    end
                catch
                    warning("Pass unrecognized attribute at ..." + ...
                        "gb(%d).CDS(%d).text, row %d",i,j,k);
                end

            elseif ~startsWith(txt(k,:),'/') % continuing prev attribute
                splt = split(txt(k,:),'"');
                attrVal{end} = [attrVal{end},splt{1}];
            end

        end

%       Write faa attributes
        faa{end+1,1} = '>';
        for k = 1:length(attr)
            faa{end,1} = [faa{end,1},'[',attr{k},'=',attrVal{k},']'];
        end

%       Write faa seq
        faa{end+1,1} = seq;

    end
end

% write faa file only when fileName is supplied
if exist('faaFileName','var')
    writecell(faa,[faaFileName,'.faa'],'QuoteStrings',false,'FileType','text');
end

end
