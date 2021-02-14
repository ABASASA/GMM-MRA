function [delta2] = ComputeGeodeticDistanceFromI_uptoscalar(A)
% A and B same szie square, semi-postaive definte, stmetric matrices
A = A / norm(A, 'fro') * sqrt(size(A,1));
egis = eig(inv(A));
% eigsInv = 1 ./ egis;

delta2 = sqrt(sum(log(egis).^2));

end

