%% Info
% This fuunction computes the Analytic first and second moments from
% the ground truth signal & distributions..
% The function is fitted to the MRA model with prjection and oulirers. 
% Inputs:
%   - rho: The ground truth shifts distibution, represented by vecotr from length L.
%   - signal: The ground truth signal, represented by vecotr from length L.
%   - projeciton: is the projection matrix which was activated on each
%                 obsevations, matrix from size Lxk.
%   - pOutlier: A scalar which represent the probablility which an
%               observation will be an outlier.
%   - CovOutlier: A matrix which represnet the outliers' covraince
%                 distibuiton.
% Output:
%   - M1: The first analytic moment.
%   - M2: The second analytic moment.
% Asaf Abas 23.08.20
function [M1, M2] = ComputeAnalyticMoments(rho, signal, sigma,...
                                projection, pOutlier, CovOutlier)

PCx = projection * circulant(signal);

% First moment
M1 = (1 - pOutlier) * PCx * rho;

% Second moment: 
M2 = (1 - pOutlier) *  PCx * diag(rho) * transpose(PCx) + ... %analytic moment of MRA
     (1 - pOutlier) * projection * sigma *  transpose(projection) + ... % bias of noise
     pOutlier * CovOutlier; % bias of outliers
end