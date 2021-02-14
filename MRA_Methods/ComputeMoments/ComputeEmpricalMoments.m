%% Info
% This fuunction computes the Emprical first and second moments from the
% observations.
% The function is fitted to the MRA model with prjection and oulirers. 
% Inputs:
%   - observations: It is a matrix which holds the observations. Each
%                   coulmn is an observation.
% Output:
%   - M1: The first empirical moment.
%   - M2: The second empirical moment.
% Asaf Abas 23.08.20
function [M1, M2] = ComputeEmpricalMoments(observations)

N = size(observations, 2); % number of observations
L = size(observations, 1); % length of each observsion (signal's length).

%% Compute the first moments
M1 = sum(observations, 2) / N;

%% Second moment
M2 = zeros(L, L);

% Compute the moment for each observaion and sum over it. 
for i = 1 : N
    M2 = M2 + observations(:, i) * transpose(observations(:, i));
end

% Avrage
M2 = M2 ./ N ; 

end
