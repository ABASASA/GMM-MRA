function [estSignal, estRho, estW, info] = GMM_L1_Iteratvie_W(startingPoints, Wstart,...
                observations, estimations, sigma, projection,pOutlier, CovOutlier,...
                numberOfWCompute, numberOfInnterIterations, FlagMeasreTime)
%% init
if ~exist('numberOfInnterIterations','var')
    numberOfInnterIterations = 2000;
elseif  numberOfInnterIterations == inf
    numberOfInnterIterations = 2000;
end

if ~exist('FlagMeasreTime', 'var')
    FlagMeasreTime = false;
end

L = size(projection, 2); % signal length
numberofStartingPoints = size(startingPoints,2);

%  data structures to save results
solutions = zeros(2 * L, numberofStartingPoints);
scores = zeros(numberofStartingPoints,1); %objective values
Infos = cell(numberofStartingPoints,1);
WsSaveVar = zeros(size(Wstart,1), numberofStartingPoints);
%% Run optimization on each guess         
for iFirstGuess = 1 : numberofStartingPoints
    currentGuess = startingPoints(:, iFirstGuess);
    currentW = Wstart;
    currentInfo = struct;
    currentInfo.CPUTime = 0;
    currentInfo.iterations = 0;
    currentInfo.ScoreVector = [inf];
    %% optimize over signal and distribution
    % The idea is to do a first computation of rho and signal. Then,
    % iteration over W -> rho & signal. I want to finish with rho and
    % signal in oreder to get it's best fit.
        
    [estSignalCurrent, estRhoCurrent, infoTmp] =...
                ComputeGMM_L1_Optimization(currentGuess,...
                currentW, observations, estimations, sigma, projection,...
                pOutlier, CovOutlier, numberOfInnterIterations, FlagMeasreTime);
            
    currentGuess = [estRhoCurrent; estSignalCurrent];
            
    % Update info
    if (FlagMeasreTime)
        currentInfo.CPUTime = currentInfo.CPUTime + infoTmp.CPUTime;
    end
    currentInfo.iterations = currentInfo.iterations + infoTmp.iterations;
    currentInfo.ScoreVector(1) = infoTmp.Score;
    for iNumWCom = 1 : numberOfWCompute
        
        %% Compute moments vector with current Est
        [momentVector] = ComputeFsMomentsQuickMRATriuOutlier(estRhoCurrent,...
                    estSignalCurrent, [], sigma,...
                    estimations, projection, pOutlier, CovOutlier);
        
        %% optimize over W
        [currentW, infoTmp] = ComputeOptimization_L1_IterativeW(momentVector,...
                            currentW, numberOfInnterIterations, FlagMeasreTime);
        
        % Update
        if (FlagMeasreTime)
            currentInfo.CPUTime = currentInfo.CPUTime + infoTmp.CPUTime;
        end
        currentInfo.iterations = currentInfo.iterations + infoTmp.iterations;
        
        
        %% optimize over signal and distribution
        [estSignalCurrent, estRhoCurrent, infoTmp] =...
                ComputeGMM_L1_Optimization(currentGuess,...
                currentW, observations, estimations, sigma, projection,...
                pOutlier, CovOutlier, numberOfInnterIterations, FlagMeasreTime);

        currentGuess = [estRhoCurrent; estSignalCurrent];

        % Update
        if (FlagMeasreTime)
            currentInfo.CPUTime = currentInfo.CPUTime + infoTmp.CPUTime;
        end
        currentInfo.iterations = currentInfo.iterations + infoTmp.iterations;
        currentInfo.ScoreVector(end+1) = infoTmp.Score;
        
    end
    WsSaveVar(:, iFirstGuess) = currentW;
    solutions(:, iFirstGuess) = currentGuess;
    Infos{iFirstGuess} = currentInfo;
    scores(iFirstGuess) = currentInfo.ScoreVector(end);
end
%% Choose the best results (via score) and return variables
[bestScore, bestScoreIndex] = min(scores);
res = solutions(:, bestScoreIndex);
info = Infos{bestScoreIndex};
info.Score = bestScore;
if ~FlagMeasreTime % to enter infinty if we didn't measure time
    info.CPUTime = inf;
end
%% Orginize output
estRho = res(1:L);
estSignal = res(L+1:end);
estW = WsSaveVar(:, bestScoreIndex);
end