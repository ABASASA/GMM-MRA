%% Info
% This function return the cost and the gradiant of GMM.
% The code itself is flexable and let to enter the moments fucntion and the
% jacobian function as an input (function handle).
% In general the GMM optimization's cost fuctnion is:
% min g(state vector, estimations)' * W * g(state vector, estimations)
% Where g is the moments function.
% Input:
% - stateVector: The state vecotor of the optimiation (in MRA it is
% [theta;signal]).
% - estimations: the estimations (to skipe the need to compute it from the
% - W: the current GMM's weights matrix.
% - sigma: the noise's coviarance matrix.
% - momentsFunction: a hangle of the moments function handle - it must reutrns the
% finel vector moments (see the function in
% GMM_Methods\Simultaneously\ComputeMomentFucntion for example).
% - jacobianFunctiton: a handle of the jacobian function for the moments
% function.
% Output:
% - cost - the results of the GMM cost function for the current state
% vector and W.
% - Grad - the current gradient of the GMM.
% 31.08.20 Asaf Abas
function [cost, Grad] = GMMCostAndGradCompute(stateVector, W,...
                                        momentsFunction, jacobianFunctiton)
%% compute cost

fSums = momentsFunction(stateVector);

cost = transpose(fSums) * W * fSums;                    

%% Compute Grad
% From https://atmos.washington.edu/~dennis/MatrixCalculus.pdf Proposition
% 12
[grad]  = jacobianFunctiton(stateVector);

Grad = transpose(fSums) * ( transpose(W) + W ) * grad;
Grad = transpose(Grad);

end