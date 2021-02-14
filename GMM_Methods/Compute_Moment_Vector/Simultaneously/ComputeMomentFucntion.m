%% Info
% This function compute the moments function (in GMM) in MRA model. 
% This function uses a pre-computed estimation of the moments (doesn't
% compute in the function). Therefore this methods is used in the
% optimization itself, not in the weights matrix computation.
%
% In general, the method compute the analytic moments with the rho and
% signal which entered to the function. Then it subtruct the moments
% estimation (emprical moments) from the analytic moments.
% Input:
% - rho: the estimated (guessed) shifts' distribution, vector from length
% L (signal's length).
% signal: the estimated signal, vector from length
% L (signal's length).
% - empricalMoment: Pre computed emprical moments.
% - projection: The projection matrix (or estiamted one).
% - pOutlier: Outliers' precent (or estiamted one), scalar.
% - CovOutlier: Outliers' covariance (or estiamted one), matrix.
% Output:
% - f - the moment vector for the current rho and signal.
% 30.08.20 Asaf Abas.
function [f] = ComputeMomentFucntion(rho, signal, sigma, empricalMoment,...
                                     projection, pOutlier, CovOutlier)
%% Get Rho and x

% Compute the analytic mometns from the current rho and signal.
[M1, M2] = ComputeAnalyticMoments(rho, signal, sigma, projection,...
                                  pOutlier, CovOutlier);
M2 = ExtractUpperTriangleMatrixVectorize(M2);% Organize M2

%% Compute the diffrence
f  = [M1 ; M2 ] - empricalMoment; 

end