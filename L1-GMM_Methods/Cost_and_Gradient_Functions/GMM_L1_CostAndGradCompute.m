%% Info
% This function return the cost and the gradiant of L1-GMM.
% The code itself is flexable and let to enter the moments fucntion and the
% jacobian function as an input (function handle).
% Input:
% - stateVector: The state vecotor of the optimiation (in MRA it is
% [theta;signal]).
% - W: the current GMM's weights matrix.
% - momentsFunction: a hangle of the moments function handle - it must reutrns the
% finel vector moments (see the function in
% GMM_Methods\Simultaneously\ComputeMomentFucntion for example).
% - jacobianFunctiton: a handle of the jacobian function for the moments
% function.
% Output:
% - cost - the results of the GMM cost function for the current state
% vector and W.
% - Grad - the current gradient of the GMM.
% 4.09.20 Asaf Abas

function [cost, Grad] = GMM_L1_CostAndGradCompute(stateVector, W,...
                                        momentsFunction, jacobianFunctiton )
%% compute cost
[fSums] = momentsFunction(stateVector);
cost = sum(abs(W .* fSums));                                 

%% Compute Grad
[grad]  = jacobianFunctiton(stateVector);

WW = repmat(W,1,size(grad,2));
Grad =  transpose(WW .* grad) * sign(W .* fSums);


end