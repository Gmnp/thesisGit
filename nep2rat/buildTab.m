load results21.mat
values = zeros([size(residss), 3]); % length(resids), min, max, mean, std
values(:,:,1) = cellfun(@length, residss);
values(:,:,2) = cellfun(@mymin, residss);
values(:,:,3) = cellfun(@mymax, residss);
% values(:,:,4) = cellfun(@mean, residss);
% values(:,:,5) = cellfun(@std, residss);

function v = mymin(z)
if isempty(z)
    v = Inf;
else
    v = min(z);
end
end

function v = mymax(z)
if isempty(z)
    v= Inf;
else
    v = max(z);
end
end



