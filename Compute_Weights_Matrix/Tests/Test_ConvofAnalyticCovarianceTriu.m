%% Info
% This test check the analytic computation of the covriance matrix. It is
% only availbe for the sceanrio without outlier, projeciton and for
% homogenius noise (represet by a scalar * Id). 
% The graph disply the decay of the error as the number of observations
% grows.
% Asaf ABas 03.09.20


%% add paths
clear;
clc;
AddPaths();

%% params
% Ground Truth input
L = 11; %Signal length

% Observations input
Ns = floor(logspace(2,5.1,6)); % Number of observations

sigmaScalar = 1; % sigma values to check
sigmaMat = sigmaScalar^2 * eye(L, L); % homogenous noise case

pOutlier = 0; % outliers precent
covOutlier = eye(L,L); % Outliers' covriance

projection = eye(L,L); % projection matrix

% Comparassion input
repTest = 12; %number of repetion for each number of observations.

% Define saving fig paramter
savingPath = 'Graphs/Test_Compare_Analytic_and_Numeric_W/';

%% Compare for each number of observations
errors = zeros(length(Ns), repTest);
for iN = 1 : length(Ns)
    N = Ns(iN); % current number of observaitons.
    tmperror = zeros(repTest,1);
    %% Define signal & distribution
    [currentGroundTruth] = GenerateStartingPoints(L, 1);
    rho = currentGroundTruth(1:L);
    signal = currentGroundTruth(L+1:end);
    [CovAn] = ComputeCovAnalyticly(signal, rho, sigmaScalar, L);

    % Repet for each numbr of observation
    parfor iRep = 1 : repTest
        
        
        %% Create observations
        observations = GenerateObservations(signal, N,...
                            rho, sigmaMat, projection, pOutlier, covOutlier);


        %% Compute Covaraince

        % Numeric
        fsFuction = @(theta, observations, sigma) ComputeMomentFucntion1By1(...
                  theta(1:L), theta((L+1) : end), observations, sigma,...
                  projection, pOutlier, covOutlier);

        [~, S] = ComputeW(currentGroundTruth, observations, sigmaMat, fsFuction);
        
        %% Compute Error
        tmperror(iRep) = norm(S(:) - CovAn(:));
    end
    errors(iN,:) = tmperror;
    disp(iN);
end
%% Display
figure;
fig = loglog(Ns, mean(errors,2), '*--');
xlabel('# Observations');
ylabel('Mean Absulote error');
title('Comparasion between: Analytic to Numeric Covriance');
% save figure
% mkdir(savingPath);
% saveas(fig, [savingPath, 'ConvarganceOfCov'],'fig'); 
% saveas(fig, [savingPath, 'ConvarganceOfCov'],'jpg'); 
