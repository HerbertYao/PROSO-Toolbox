% Parse HumanCyc .col data struct into a table
function dat_out = datParse(dat_in,subjects,features)

if ~exist('features','var')
    features = cell(length(subjects),1);
end

dat_out = cell(1,length(subjects));
idx = 1;

for i = 1:length(dat_in)
    if strcmp(dat_in{i,1},'//')
        idx = idx + 1;
        dat_out{idx,1} = '';
    elseif any(strcmp(subjects,dat_in{i,1}))
        col = find(strcmp(subjects,dat_in{i,1}));

        for j = 1:length(col)
            if isempty(features{col(j)})
                if isempty(dat_out{idx,col(j)})
                    dat_out{idx,col(j)} = dat_in{i,2};
                else
                    dat_out{idx,col(j)} = [dat_out{idx,col(j)},',',dat_in{i,2}];
                end
                break;

            elseif contains(dat_in{i,2},features{col(j)})
                if isempty(dat_out{idx,col(j)})
                    dat_out{idx,col(j)} = dat_in{i,2};
                else
                    dat_out{idx,col(j)} = [dat_out{idx,col(j)},',',dat_in{i,2}];
                end
                break;
            end
        end
    end

end

end