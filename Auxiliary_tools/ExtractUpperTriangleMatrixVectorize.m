%% Info
% This function extrats the upper trianglr part from an matrix and
% vectorize it.
% In general using triu() function returns the matrix with 0's in non upper
% triangle's entries. However, for our use (computing W) we want to the
% lower left. Therefore we extract the indeces of the upper matrix
% (including the diagonal) and take the relevent entries in the matrix.
% Input:
%   - A: A square matrix.
% Output:
%   - ATUp: The vectorize upper triangle part of A.
% Asaf Abas 23.08.20

function [ATUp] = ExtractUpperTriangleMatrixVectorize(A)
n = size(A,1);
onesMat = true(n,n);
onesMat = triu(onesMat);
ATUp = A(onesMat);
end