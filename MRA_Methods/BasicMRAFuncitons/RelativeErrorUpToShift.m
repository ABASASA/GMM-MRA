%% Info
% This funciton compare between 2 vector (L2- wise) up to
% discreate shift.
% The code basiclly try all combinations of shift of one of the
% vector(XEst) and compare it to the other (x). Then, choose the shift with
% the smallest L2 error and return it.
% input:
% - x and xEst to vector from the same size which we will want to find the
% smallest rellative error between them.
% Output:
% - error - the smallest rellative error, L2 - wise and up to shift.
% 30.8.20 Asaf Abas
function [error] = RelativeErrorUpToShift(x, xEst)
%% Check the sizes
if (size(xEst,1) ~= 1)
    xEst = transpose(xEst);
end
if (size(x,1) ~= 1)
    x = transpose(x);
end
%% Compare
% create all shifts of xEst
circulentxEst = circulant(xEst); % Make sure is a row
X = repmat(x, length(x),1); % Copies of x
%% find minimumm
error = min(real(sqrt(sum((circulentxEst - X) .^ 2, 2 )) / norm(x,2)));

end