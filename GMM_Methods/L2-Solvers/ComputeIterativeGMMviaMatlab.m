%% Info
% iterative GMM, after each iteration (except the first) updates the
% weights matrix.
% The opimization itself is computed by another fuctnion, this is just a
% function to mange the organization of the iteration and first guesses.
% It chooses the estimation with the lowest score (after all the step
% for each first guess).
% This is support: prjection and outlieres.
% Input:
% startingPoints - a matrix where each column is a first guess ([rho;signal])
% WFirst - the first iteration's weights matrix.
% observations - needed to compute the weights matrix.
% empricalMoment - emprical moments.
% sigma - noise's covriance.
% momentFuction1by1 - function handle which needed to compute weight
% matrix.
% projection - projection matrix
% pOutlier - precent of outliers.
% CovOutlier - covriance of the outliers.
% maxNumOfIterations - max number of opimization iterations.
% numberOfGMMIterations - number of GMM iterations.
% FlagMeasreTime - true if you want to save cou time.
% Output:
% estSignal - estimated signal (best score).
% estRho - estimated shift distibution (best score).
% ASaf Abas 31.08.20

function [estSignal, estRho, info, moreData] = ComputeIterativeGMMviaMatlab(startingPoints,...
                WFirst, observations, empricalMoment, sigma,...
                momentFuction1by1, projection, pOutlier, CovOutlier,...
                maxNumOfIterations, numberOfGMMIterations, FlagMeasreTime)
%% init
L = size(startingPoints, 1) ./ 2; % signals length
numberofStartingPoints = size(startingPoints,2); % number of first guess

% data matrix to sav results in each iteration for each guess
solutionsSignal = zeros(L, numberofStartingPoints, numberOfGMMIterations);
solutionsRho = zeros(L, numberofStartingPoints, numberOfGMMIterations);

scores = zeros(numberofStartingPoints, numberOfGMMIterations); % cost function scores.
Infos = cell(numberofStartingPoints, 1);
%% Optimization
%  run over all first guesses
for iRep = 1 : numberofStartingPoints
    
    firstGuess = startingPoints(:, iRep); % choose current first guess
    
    numberOfIterations = 0; 
    if FlagMeasreTime
        CPUTime = 0;
    end
    %% First Step
    [estSignal, estRho, currentInfo] = OptimizationIterationGMMinMRAviaMatlab(firstGuess,...
            WFirst, empricalMoment, sigma, projection, pOutlier, CovOutlier,...
                        maxNumOfIterations, FlagMeasreTime);
                    
    currentGuess = [estRho; estSignal];
    solutionsRho(:, iRep, 1) = estRho;
    solutionsSignal(:, iRep, 1) = estSignal;
    numberOfIterations = currentInfo.iterations;
    if FlagMeasreTime
        CPUTime = currentInfo.CPUTime;
    end
    scores(iRep, 1) = currentInfo.Score;
    %% Iterative Steps
    for iStep = 2 : numberOfGMMIterations
        
        % Update Weights matrix
        if FlagMeasreTime
            tic
        end
       [currntW, ~] = ComputeW(currentGuess , observations, sigma, momentFuction1by1);
        if FlagMeasreTime
            timeW = toc;
            CPUTime = timeW + CPUTime;
        end
        
        % Compute next step
        [estSignal ,  estRho, currentInfo] = ...
            OptimizationIterationGMMinMRAviaMatlab(currentGuess, currntW, ...
                     empricalMoment, sigma, projection, pOutlier, CovOutlier,...
                     maxNumOfIterations, FlagMeasreTime); 
         if FlagMeasreTime
            CPUTime = currentInfo.CPUTime + CPUTime;
        end
        
        currentGuess = [estRho; estSignal];
        %% Save Data
        solutionsSignal(:,iRep, iStep) = estSignal;
        solutionsRho(:,iRep, iStep) = estRho;
        scores(iRep, iStep) = currentInfo.Score;
        numberOfIterations = currentInfo.iterations + numberOfIterations;
    end
    currentInfo.iterations = numberOfIterations;
    currentInfo.CPUTime = CPUTime;
    Infos{iRep} = currentInfo;

end
%% Find Best Score
[~,bestScoreIndex] = min(squeeze(scores(:,end)));
info = Infos{bestScoreIndex};
%% Orginize output
estRho = solutionsRho(:, bestScoreIndex, end);
estSignal = solutionsSignal(:, bestScoreIndex, end);

moreData  = struct;
moreData.scores = squeeze(scores(bestScoreIndex,:));
moreData.rho = squeeze(solutionsRho(:,bestScoreIndex,:));
moreData.singal = squeeze(solutionsSignal(:,bestScoreIndex,:));
end