%% Info
% The function preform iterative-L1-GMM with optimization over weights
% vector (in L1 it is a vector). 
% startingPoints - a matrix where each column is a first guess ([rho;signal])
% WFirst - the first iteration's weights matrix.
% observations - needed to compute the weights matrix.
% empricalMoment - emprical moments.
% sigma - noise's covriance.
% projection - projection matrix
% pOutlier - precent of outliers.
% CovOutlier - covriance of the outliers.
% maxNumOfIterations - max number of opimization iterations.
% numberOfGMMIterations - number of GMM iterations.
% FlagMeasreTime - true if you want to save cou time.
% Output:
% estSignal - estimated signal (best score).
% estRho - estimated shift distibution (best score).

function [estSignal, estRho, info, moreData] = Iterative_GMM_L1_Optimization(...
        startingPoints, WFirst, empricalMoment, sigma,...
        projection, pOutlierEstimator, CovOutlier, maxNumOfIterations,...
        numberOfGMMIterations, FlagMeasreTime, groundTruth)
%% Init
numberofStartingPoints = size(startingPoints,2);
L = size(startingPoints, 1) ./ 2; % signal's length

solutionsSignal = zeros(L, numberofStartingPoints, numberOfGMMIterations);
solutionsRho = zeros(L, numberofStartingPoints, numberOfGMMIterations);

scores = zeros(numberofStartingPoints, numberOfGMMIterations); %obective values
Infos = cell(numberofStartingPoints, 1);
%% Going over all initial guesses
for iRep = 1 : numberofStartingPoints
    
    firstGuess = startingPoints(:, iRep);
    numberOfIterations = 0;
    if FlagMeasreTime
        CPUTime = 0;
    end
    
    
    %% First Step
    [estSignal, estRho, currentInfo] = OptimizationIteration_L1_GMM_viaMatlab(...
            firstGuess, WFirst, empricalMoment, sigma,...
                projection, pOutlierEstimator, CovOutlier,...
                maxNumOfIterations, FlagMeasreTime);
    nextGuess = [estRho; estSignal];
    solutionsRho(:, iRep, 1) = estRho;
    solutionsSignal(:, iRep, 1) = estSignal;
    numberOfIterations = currentInfo.iterations;
    scores(iRep, 1) = currentInfo.Score;
    if FlagMeasreTime
        CPUTime = CPUTime + currentInfo.CPUTime;
    end
    
    
    Wold = WFirst;
    %% Iterative Steps
    for iStep = 2 : numberOfGMMIterations
        if FlagMeasreTime
        	tic;
        end
        % Compute the moment vector for current estimation.

        if exist('groundTruth', 'var')
            
            momentVector = ComputeMomentFucntion(groundTruth(1:L),...
                         groundTruth(L+1:end), sigma, empricalMoment,...
                         projection, pOutlierEstimator, CovOutlier);
        else
            momentVector = ComputeMomentFucntion(estRho,...
                         estSignal, sigma, empricalMoment,...
                         projection, pOutlierEstimator, CovOutlier);
        end
        % Optimize weights vecotr
       [Wcurrent, ~] = Optimization_L1_Find_W_viaManopt(momentVector,...
                  Wold, maxNumOfIterations, FlagMeasreTime);
        if FlagMeasreTime
            timeW = toc;
            CPUTime = timeW + CPUTime;
        end
        % Compute next guess
        [estSignal, estRho, currentInfo] = OptimizationIteration_L1_GMM_viaMatlab(...
            nextGuess, Wcurrent, empricalMoment, sigma,...
                projection, pOutlierEstimator, CovOutlier,...
                maxNumOfIterations, FlagMeasreTime);
        if FlagMeasreTime
            CPUTime = CPUTime + currentInfo.CPUTime;
        end
        nextGuess = [estRho; estSignal];
        %% Save Data
        solutionsSignal(:,iRep, iStep) = estSignal;
        solutionsRho(:,iRep, iStep) = estRho;
        scores(iRep, iStep) = currentInfo.Score;
        numberOfIterations = currentInfo.iterations + numberOfIterations;
        
        Wold = Wcurrent;

    end
    currentInfo.iterations = numberOfIterations;
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