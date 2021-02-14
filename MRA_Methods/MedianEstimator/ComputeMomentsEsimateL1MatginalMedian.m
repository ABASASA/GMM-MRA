function [M1, M2] = ComputeMomentsEsimateL1MatginalMedian(observations, sigma, projection,divideBy)

N = size(observations, 2);
L = size(observations, 1);

M1 = median(observations, 2);
M2s = zeros(L^2, N);
M2 = zeros(L, L);
%% Compute for each sample 
for i = 1 : N
    a = observations(:,i) * transpose(observations(:, i));
    
    M2s(:,i) = a(:);
end
%% Add projection and un-bias
M2 = median(M2s,2);
M2 = reshape(M2,[L,L]) -  projection * sigma *  transpose(projection) ; 

end
