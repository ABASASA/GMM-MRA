%% Info
% This function compute the variance of the first moment analyticly.
% Asaf Abas 4.9.20
function [cov] = ComputeVarFirstMomentsAnaylticly(signal, rho, sigma, L)

%% Compute the moments
[M1, M2] = ComputeAnalyticMoments(rho, signal, sigma.^2 * eye(L,L),...
                                eye(L,L), 0, eye(L,L));

%% Compute First Moment's cov
cov = (M2  - M1 * (M1'));

end