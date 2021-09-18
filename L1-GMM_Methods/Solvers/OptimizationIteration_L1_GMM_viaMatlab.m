%% Info
% This function preforms a L1-GMM optimization using matlab's contranted
% optimizer.
% The function defines the constraint and setting of the optimizatiton.
% Also get the weight mattrix and the intial guess (starting guess) and
% define the moment fcunton and jacobian which the optimization will use.
% Input:
% startingPoint - the intial guess
% W - the  weights matrix.
% empricalMoment - emprical moments.
% sigma - noise's covriance.
% projection - projection matrix
% pOutlier - precent of outliers.
% CovOutlier - covriance of the outliers.
% maxNumOfIterations - max number of opimization iterations.
% FlagMeasreTime - true if you want to save cou time.
% Output:
% estSignal - estimated signal
% estRho - estimated shift distibution
% ASaf Abas 04.09.20
function [estSignal, estRho, info] = OptimizationIteration_L1_GMM_viaMatlab(...
                startingPoint, W, empricalMoment, sigma,...
                projection, pOutlier, CovOutlier,...
                maxNumOfIterations, FlagMeasreTime)

if ~exist('maxNumOfIterations','var')
    maxNumOfIterations = 2000;
elseif  maxNumOfIterations == inf
    maxNumOfIterations = 2000;
end

if ~exist('FlagMeasreTime', 'var')
    FlagMeasreTime = false;
end

L = length(startingPoint) / 2;
%% Enters the Cost & Grad's Functions
momentsFunction = @(theta)...
                         ComputeMomentFucntion(theta(1:L),...
                         theta((L+1) : end), sigma, empricalMoment,...
                         projection, pOutlier, CovOutlier);              

jacobianFunctiton = @(theta) ComputeJacobianMRA(theta(1:L),...
                                            theta((L+1) : end), projection, pOutlier);

%% Constraints - without Limits on rho
A = [-eye(L), zeros(L,L)];
b = [zeros(L,1)];
Aeq =[ones(1,L), zeros(1,L)];
beq = 1;

%% Solver Options - Copied from Nir's code
options = optimoptions(@fmincon,'Algorithm','interior-point');
options.MaxIterations          = maxNumOfIterations;
options.MaxFunctionEvaluations = 10000;
options.FunctionTolerance = 1e-25;
options.ConstraintTolerance =  1e-10;
options.StepTolerance =  1e-16; 
options.Display = 'off';
options.SpecifyObjectiveGradient = true(1);
options.OptimalityTolerance = 1e-20; %NIR
options.FiniteDifferenceStepSize = 1e-25;

%% Define the optimization Function

optimizationFunction = @(x) GMM_L1_CostAndGradCompute(x, W,...
                                momentsFunction, jacobianFunctiton);

%% Optimization
if FlagMeasreTime
    tic;
end

% Optimization
[solution, score, flag1, info] = fmincon(optimizationFunction,...
                startingPoint, A, b, Aeq, beq, [], [], [], options);

if FlagMeasreTime
    CPUTime = toc; 
end

info.Score = score;
info.Flag = flag1;
if FlagMeasreTime
    info.CPUTime = CPUTime;
else
    info.CPUTime = inf;
end
%% Orginize output
estRho = solution(1 : L);
estSignal = solution(L+1 : end);                            
end