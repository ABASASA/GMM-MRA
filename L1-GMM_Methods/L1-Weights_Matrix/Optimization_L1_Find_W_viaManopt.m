%% Info
% This function preforms an optimization in order to find optimal weights
% vector for L1 cost function for a known moment vector.
% The fucntion uses Manopt to solve: argmin_w (norm(W .* moment vector,1))
% Inputs:
% momnetVector - the moment vector - fixed durring the optimization.
% Wstart - the intialguess for weights vector.
% maxNumOfIterations - max number of opimization iterations.
% FlagMeasreTime - true if you want to save cou time.
% output:
% finalW - the opimization's results.
% info -  struct to hold data as score and time.
function [finalW, info] = Optimization_L1_Find_W_viaManopt(momnetVector,...
                            Wstart, maxNumOfIterations, FlagMeasreTime)

if ~exist('maxNumOfIterations','var')
    maxNumOfIterations = 2000;
elseif  maxNumOfIterations == inf
    maxNumOfIterations = 2000;
end

if ~exist('FlagMeasreTime', 'var')
    FlagMeasreTime = false;
end

Wsize = length(Wstart);

%% Define fucntions
costFunc = @(W) sum(abs(W .* momnetVector));
gradFunc = @(W) sign(W .* momnetVector) .* momnetVector;

%% Define manifold's problem
manifold = spherefactory(Wsize, 1);
problem.M = manifold;

problem.cost = costFunc;
problem.egrad = gradFunc;
%% Options
options.maxiter = maxNumOfIterations;
options.verbosity = 0;
%% Optimization
if FlagMeasreTime
        tic;
end
% Optimization
% checkgradient(problem);
warning('off', 'manopt:getHessian:approx')
[finalW, score, infoManopt, ~] = rlbfgs(problem, Wstart, options);

% optional measutre time
if FlagMeasreTime
    CPUTime = toc; 
end
% save info
info = struct;
info.iterations =  length([infoManopt.iter]);
info.Score = score;
if FlagMeasreTime
    info.CPUTime = CPUTime;
end

end
