%% Info
% This code generates observations to MRA + projection + outliers.
% The basic idea is first to create the non-outliers, AKA ok observations.
% Each OK observation is created following the next procidure:
% 1. Sample a shift of the original signal according to rho (shift's
% ditibution).
% 2. Add noise according to sigma (with Normal(0, sigma)) - sigma is a
% matrix from size (LxL where L is the signal's length).
% 3. Project the observation with the projection matrix (input:
% "projection"). The projection matrix is from size (LxK - where K is thhe
% observations' size after the projaction - usally smaller).
% 
% Outliers created according an simpler procidure:
% 1. Gussian noise accoding to sigmaOutlier matrix from size( LxK).
%
% We decies if each obervation is OK or an outliers according to pOutlier (a scalar).
% For each observation sampling a uniform dist(Uni([0,1})) if it bigger then pOutlier
% it is an ok observation, otherwise an outlier.
% Input:
% - signal - The ground truth signal, a vector from length L.
% - numberOfObservations - The total number of observation we want to
% create (a scalar).
% - rho - The shift's distibution, a vector from length L.
% - sigma - The covriance matrix of the noises' gaussan distibution of the
% ok observations. The matrix is from size LxL.
% - projaction - The projaction matrix, from size LxK. In order not to have
% projection enter eye(L,L).
% - pOutlier - The precent of outliers in the code (a scalar). In order not
% to have outlier enter 0.
% - sigmaOutliers - the std of the guassain distibuiton of the
% outliers.
% Output:
% - observations - The matrix will the observations (in each columns). size
% [K x numberOfObservations].
% 30.08.20 Asaf Abas
function [observations] = GenerateObservations(signal, numberOfObservations,...
                rho, sigma, projection, pOutlier, sigmaOutlier)


L = length(signal); % Signal's length

observations = zeros(L, numberOfObservations); % initialize variable

CDFRho = cumsum(rho); % compute the CFD of the distribuiton

Csignal = circulant(signal); % circulant matrix of the signal
%% Compute which entries will be ouliers
% For each observation I decide if it an outlier or an regular observation.
% The process is to for each observation to sample from Uni([0,1]) and
% check if it is samller than the precent of outliers (pOutlier) in the scenairo.
% Then make lists of the indeces of observation with are OK and outliers

% Randomize if obsevration is an outleir.
isOulier = rand(numberOfObservations, 1) < pOutlier; 

indecesRegulearMRA = find(~isOulier); % OK observations list
indecesOutlier = find(isOulier); % Outliers observations list
%% Generatet OK obsevations
% For each OK observation we first shift the signal (accrding to the
% shifts' distibution - rho). Than, adding a guassion noise, then
% projection the data using the projection matrix.

% Choose the shift for each signal. The shift is sampled by sampling a
% variable from Uni([0,1]) then cheking which shift it's cronspond through
% the invarse of the CFD of rho.
for i = 1 : length(indecesRegulearMRA)
    tmp = rand;
    jShift = find(tmp <= CDFRho, 1);
    observations(:, indecesRegulearMRA(i) ) = Csignal(:, jShift);
end

% Adding noise - guassion according to the covraince matrix sigma.
observations(:,indecesRegulearMRA) = observations(:,indecesRegulearMRA) + ...
        transpose(mvnrnd(zeros(L,1), sigma, length(indecesRegulearMRA)));

% Project the observations
observations = projection * observations;

%% Generate outliers (override noise + projection)

observations(:, indecesOutlier ) = sigmaOutlier *...
            randn(size(projection,1), length(indecesOutlier));
end