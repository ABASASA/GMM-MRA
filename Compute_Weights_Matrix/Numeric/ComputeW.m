%% Info
% Compute the weights matrix for GMM for some state vector.
% In the case of the MRA the optimal wieght is the invarse of the covriance
% matrix.
% The method compute the moment funciton in for every observation them
% compute the covriance matrix.
% Input:
% - stateVector: the state vector which will for it computed the analytic
% moments ([rho;signal]) (2L length vector).
% - observations: the matrix of observaitons.
% - sigma: the noise's covirance matrix.
% - momentFuction1by1: a handle to the moment function which take each time
% 1 observations only.
% Output:
% - W: the weight matrix, the covriance's invarse.
% - covMatrix: the covriance matrix.
% Asaf Abas 31.08.20
function [W, covMatrix] = ComputeW(stateVector, observations, sigma, momentFuction1by1)

% Compute the moment vector for each observation 
[fs] = momentFuction1by1(stateVector, observations, sigma);

% compute covairance
covMatrix = cov(transpose(fs));

% compue invaes
W = pinv(covMatrix);

end