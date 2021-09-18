%% Info
% This function compute the moments function (in GMM) in MRA model by 1
% by 1. 1-by-1 means that e compute each moments (for each observation)
% sepertly and not return it's mean. This is used for computing the
% covraicne matrix (and therrefore weight matrix).

% In general, the method compute the analytic moments with the rho and
% signal which entered to the function. Then it subtruct the moments
% estimation (emprical moments) of each observation sepertly.
% Input:
% - rho: the estimated (guessed) shifts' distribution, vector from length
% L (signal's length).
% signal: the estimated signal, vector from length
% L (signal's length).
% - sigma: the covriances matrix of the noise
% - projection: The projection matrix (or estiamted one).
% - pOutlier: Outliers' precent (or estiamted one), scalar.
% - sigmaOutlier: Outliers' covariance (or estiamted one), matrix.
% Output:
% - fs - this is a matrix, where each coulmn is the moment vector for each
% observations.
% 3.09.2.0 Asaf Abas.
function [fs] = ComputeMomentFucntion1By1DirectMoMent3(rho, signal, observations,...
                sigma, projection, pOutlier, sigmaOutlier)
%% General params
L = size(observations, 1);
N = size(observations, 2);
[indecesM3] = ChooceIndecesM3(L);

%% Compute analytic moments
[M1, M2] = ComputeAnalyticMoments(rho, signal, sigma, projection, pOutlier, sigmaOutlier);
M2 = ExtractUpperTriangleMatrixVectorize(M2);

%% Intial Saveing data object
fOutputLength = length(M1) + length(M2) + size(indecesM3,1);
fs = zeros(fOutputLength, N);

%% Compute the moment vector for each observations seperly
for i = 1 : N
    
    [M1Est, M2Est] = ComputeEmpricalMoments(observations(:,i));
    M2Est = ExtractUpperTriangleMatrixVectorize(M2Est);
    [M3Emp] = ComuteM3Empric(observations(:,i),indecesM3);

    fs(:,i) = [(M1Est); (M2Est) ; M3Emp];
end
end