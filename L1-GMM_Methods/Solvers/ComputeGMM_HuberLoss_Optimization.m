function [estSignal, estRho, info] = ComputeGMM_PseudoHuberLoss_Optimization(startingPoints, W,...
                observations, estimations, sigma, projection,pOutlier, CovOutlier,...
                numberOfLS2GMMIter, FlagMeasreTime)
global itNum;
itNum = 1;
if ~exist('numberOfLS2GMMIter','var')
    numberOfLS2GMMIter = 2000;
elseif  numberOfLS2GMMIter == inf
    numberOfLS2GMMIter = 2000;
end

if ~exist('FlagMeasreTime', 'var')
    FlagMeasreTime = false;
end
L = size(projection, 2);
numberofStartingPoints = size(startingPoints,2);

%% Enters the Cost & Grad's Functions
fsFuctionQuick = @(theta, observations, sigma, estimations)...
                     ComputeFsMomentsQuickMRATriuOutlier(theta(1:L),...
                         theta((L+1) : end), observations, sigma,...
                         estimations, projection, pOutlier, CovOutlier);                

gradFsFunction = @(theta) (1 - pOutlier) * ComputeGradMRATriu(theta(1:L),...
                                            theta((L+1) : end), projection);

%% Constraints - without Limits on rho
A = [-eye(L), zeros(L,L)];
b = [zeros(L,1)];
Aeq =[ones(1,L), zeros(1,L)];
beq = 1;

%% Solver Options - Copied from Nir's code
options = optimoptions(@fmincon,'Algorithm','interior-point');
options.MaxIterations          = numberOfLS2GMMIter;
options.MaxFunctionEvaluations = 10000;
options.FunctionTolerance = 1e-25;
options.ConstraintTolerance =  1e-10;
options.StepTolerance =  1e-16; 
options.Display = 'off';
options.SpecifyObjectiveGradient = true(1);
options.OptimalityTolerance = 1e-20; %NIR
options.FiniteDifferenceStepSize = 1e-25;
% options.PlotFcn = [];% @optimplotfval ;
% options.CheckGradients = true;
%% Define the optimization Function

optimizationFunction = @(x) GMM_PseudoHuberLoss_CostAndGradCompute(x, observations, W,...
                                sigma, fsFuctionQuick, gradFsFunction, estimations);
%% 
solutions = zeros(2 * L, numberofStartingPoints);
scores = zeros(numberofStartingPoints,1); %objective values
Infos = cell(numberofStartingPoints,1);
CPUTime = zeros(numberofStartingPoints,1);

%% Optimization
for iRep = 1 : numberofStartingPoints
    % optional measutre time
    if FlagMeasreTime
        tic;
    end
    % Optimization
    [solutions(:,iRep), scores(iRep), ~, Infos{iRep}] = fmincon(optimizationFunction,...
        startingPoints(:, iRep), A,b,Aeq,beq,[],[],[],options);
    % optional measutre time
    if FlagMeasreTime
        CPUTime(iRep) = toc; 
    end
end

[bestScore, bestScoreIndex] = min(scores);
res = solutions(:, bestScoreIndex);
info = Infos{bestScoreIndex};
info.Score = bestScore;
if FlagMeasreTime
    info.CPUTime = CPUTime(bestScoreIndex);
else
    info.CPUTime = inf;
end
%% Orginize output
estRho = res(1:L);
estSignal = res(L+1:end);

% 
% figCond1 = figure(figCond1)
% hold on
% semilogy(cond1, '*-');
% 
% figCond2 = figure(figCond2)
% hold on
% semilogy(cond2, '*-');
% figure;
% semilogy(history,'r*-');
end