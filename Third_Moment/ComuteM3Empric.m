function [M3] = ComuteM3Empric(observations, indeces)


N = size(observations, 2); % number of observations
L = size(observations, 1); % length of each observsion (signal's length).

numInd = size(indeces,1);
M3 = zeros(numInd , 1);

for ind = 1 : numInd
    k1  = indeces(ind, 1);
    k2  = indeces(ind, 2);
    k3  = indeces(ind, 3);

    xtmp1 = observations(k1, :);
    xtmp2 = observations(k2, :);
    xtmp3 = observations(k3, :);
    M3(ind) = sum(xtmp1 .* xtmp2 .* xtmp3);

end

M3 = M3 ./ N;
end