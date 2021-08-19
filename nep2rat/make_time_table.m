load results21.mat
nprobs = size(values,1);
nbalgs = 5;

fileID = fopen('tabTimings.txt','w');
[~,minTim] = min(timings,[],2);
[~,maxTim] = max(timings,[],2);

for i=1:nprobs
    nep{i} = replace(nep{i}, '_', '\_');
    fprintf(fileID, '\\texttt{%s (%d)} &', nep{i}, nevs(i));
    for k=1:nbalgs
        if timings(i,k) > 10
            fprintf(fileID, '\\rc %7.2e', timings(i,k));
        else
            fprintf(fileID, '%7.2e', timings(i,k));
        end
        if k < nbalgs
            fprintf(fileID,' &');
        end
    end

    fprintf(fileID, '\\\\\n');
end

fclose(fileID);