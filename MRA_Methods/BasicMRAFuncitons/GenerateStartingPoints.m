%% Info
% This function create a strarting points of the optimization and the grund
% truth. 
% Each points is contructed from 2 parts: distibuiton and signal.

% The distiubtion is a vector of length L, which each entry is sampled
% from uniform distibution between [0,1]. After sampling the entry we
% devide each of them by the total sum of the entries. This way the total
% sum is 1.

% The signal is a vector of length L. Each entry is smapled iid from
% guassion distribution and then the vector is normlize.

% Input:
%   - L: The signal's length. A scalar.
%   - numberOfPoints: number of points we want to generate.
%
% Output:
%   - startingPoints: A matrix where each coulmn is a Starting point (or
%   ground truth) [distribuiton; signal].
% Asaf Abas 23.08.20
function [startingPoints] = GenerateStartingPoints(L, numberOfPoints)

% This matrix will hold the generated points.
startingPoints = zeros(2 * L, numberOfPoints); 


for i = 1 : numberOfPoints % For each point
    % Distubition
    rho = rand(L,1); 
    rho = rho ./ sum(rho);
    % Signal
    signal = randn(L,1);
    signal = signal / norm(signal);

    startingPoints(:, i) = [rho; signal]; % Save
end
end