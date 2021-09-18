%% add paths
clear;
clc;
AddPaths();
%% params
% Ground Truth input
L = 10; %Signal length

% Observations input
numberOfObservations = 100000; % Number of observations
sigmaScalar = 0.5 ;



sigmaDiag = ones(L, 1); % homogenous noise case
sigma = (sigmaScalar^2) * diag(sigmaDiag);

pOutlier = 0; % outliers precent
sigmaOutlier = 10 * eye(L,L);
projection = eye(L,L);


% Comparassion input
numberOfStartingPoints = 1; % number of intial guesses

% Optimization input
maxNumOfIterations = inf; % number of iteration in each GMM optimization
numberOfStepsGMM = 1; % number of steps in the GMM

covOurlier = sigmaOutlier .^ 2;
    
    
%% Sample a ground truth (signal + shifts' distibution)

[currentGroundTruth] = GenerateStartingPoints(L, 1);

rho = currentGroundTruth(1:L);
signal = currentGroundTruth(L + 1 : end);



%% Create observations
observations = GenerateObservations(signal, numberOfObservations,...
                    rho, sigma, projection, pOutlier, sigmaOutlier);

[indecesM3] = ChooceIndecesM3(L);

[M1Est, M2Est] = ComputeEmpricalMoments(observations);
        M2Est = ExtractUpperTriangleMatrixVectorize(M2Est);
[M3Est] = ComuteM3Empric(observations,indecesM3);
empricalMoment = [M1Est; M2Est;M3Est];

%% Random a starting point for the optimization
[startingPoints] = GenerateStartingPoints(L, numberOfStartingPoints);

%% Define computation method to W
momentFuction1by1 = @(theta, observations, sigma) ComputeMomentFucntion1By1(...
          theta(1:L), theta((L+1) : end), observations, sigma,...
          projection, pOutlier, covOurlier);

momentFuction1by1Direct = @(theta, observations, sigma) ComputeMomentFucntion1By1DirectMoMent3(...
                  theta(1:L), theta((L+1) : end), observations, sigma,...
                  projection, pOutlier, covOurlier);   
%% Compute W

W = eye(length(empricalMoment), length(empricalMoment));

WGMM = ComputeW(currentGroundTruth , observations, sigma, momentFuction1by1Direct);


%% Enters the Cost & Grad's Functions
fsFuctionQuick = @(theta, sigma1, estimations) ComputeMoment3Fucntion(...
     theta.Rho, theta.Signal, sigma1, empricalMoment, projection,...
        pOutlier, covOurlier, indecesM3, sigmaScalar);       

gradFsFunction = @(theta) ComputeJacobianM3MRA(theta.Rho, theta.Signal,...
                        projection, pOutlier, sigmaScalar);


optimizationFunction = @(x) GMMCostAndGradComputeManopt(x, empricalMoment, W,...
                                sigma, fsFuctionQuick, gradFsFunction);
costFun = @(x) nth_output(1,optimizationFunction,x);                            
egradFun = @(x) nth_output(2,optimizationFunction,x); 

%% Define the manifold - to manopt
elements = struct();
elements.Signal = euclideanfactory(L);
elements.Rho = multinomialfactory(L);
manifold = productmanifold(elements);
problem.M = manifold;

% Input cost & gradient funciton
problem.cost = costFun;
problem.egrad = egradFun;
%% options
options.maxiter = 1000;
options.verbosity = 0;
options.tolgradnorm = 1e-10;
options.minstepsize = 1e-11;
options.totalcost = 1e-10;

%% Optimization

tmpStartingPoint.Rho =  startingPoints(1:L);
tmpStartingPoint.Signal = startingPoints(L + 1: end);

checkgradient(problem);
[tmpSol, scores, Infos, ~] = steepestdescent(problem, tmpStartingPoint, options);
%%
figure;
subplot(2,1,1)
plot(tmpSol.Signal);
subplot(2,1,2);
plot(signal)
[tmpSol.Signal, signal]

