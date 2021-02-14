%% Info
% This function compute analytilcy the covarance matrix of the moment function.
% Notice tis function is for the case when the sceond moment is only he
% upper half triangle.
% It works only for non-outliers and non-projection cases  with homogenous
% noise (represented by a scalar).
% The function compute each part speprtly: Var(M1), cov(M1, M2) and
% Var(M2).
% Asaf Abas 3.09.20

function [CovAn] = ComputeCovAnalyticly(signal, rho, sigmaScalar, L)

varAnM2 = ComputeVarSecondMomentsAnayliticly(signal, rho, sigmaScalar, L);
covAnM2M1 = Compute_Cov_First_Second_MomentsAnayliticly(signal, rho, sigmaScalar, L);
varAnM1 = ComputeVarFirstMomentsAnaylticly(signal, rho, sigmaScalar, L);

CovAn = [varAnM1, transpose(covAnM2M1);
         covAnM2M1, varAnM2];
%  CovAn = varAnM2;
end