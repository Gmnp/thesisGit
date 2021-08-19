buildTab;
nprobs = size(values,1);
nbalgs = 5;

fileID = fopen('tabExp.txt','w');

for i=1:nprobs
    nep{i} = replace(nep{i}, '_', '\_');
    fprintf(fileID, '\\texttt{%s (%d)} &', nep{i}, nevs(i));
    for k=1:nbalgs
        for j = 1:3
            if j==1 % we have integers       
                if values(i,k,j) == nevs(i)
                    fprintf(fileID, '\\multicolumn{1}{c}{%d}', values(i,k,j));
                else
                    fprintf(fileID, '\\multicolumn{1}{c}{\\rc %d}', values(i,k,j));
                end
            else
                fprintf(fileID, '\\multicolumn{1}{c}{%7.0e}', values(i,k,j));
            end
            if k < nbalgs  || j < 3
                fprintf(fileID,' &');
            end
        end
    end
    fprintf(fileID, '\\\\\n');
end

fclose(fileID);