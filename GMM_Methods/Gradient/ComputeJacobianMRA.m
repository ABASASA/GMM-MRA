%% Info
% Computes the Jacobian matrix for MRA model with:
% - projection matrix
% - only upper tirangle of the second moment
% - prenect of outliers
% 31.08.20 Asaf Abas
function [jacobian] = ComputeJacobianMRA(rho, signal, projection, pOutlier)

L = length(signal); % Origianl signal's length
k = size(projection, 1); % observation's length after projection.

%% pre-compute matrices to reduce computations
Crho = circulant(rho);
Cx = circulant(signal);

PCx = projection * Cx;
CxDrho = Cx * diag(rho);

%% Create indces for the M2 (after vectorization)
k1s = repmat([1 : k], 1, k);
k2s = repmat([1:k], k, 1);
k2s = k2s(:);
%% Compute first moment

jacobianM1rho = PCx;
jacobianM1X = projection * Crho;
jacobianM1 =  [jacobianM1rho,  jacobianM1X];

%% Compute sescond moment

% Define variables
jacobianM2rhoTmp = zeros((k + 1) * k / 2, L); 
jacobianM2XTmp = zeros((k + 1) * k / 2, L); 

% We take only the upper triangle of te second moment and
% vectorize it. Therefore we need to track which index are we (in the
% vectorze seocnd moment - "currentIndex").
currentIndex = 0; 
for index = 1 : length(k1s)

    k1 = k1s(index); % compute current index in matrices indeces
    k2 = k2s(index); % compute current index in matrices indeces
    
    if(k1 > k2) % Check if it is in triu
        % the indeces in the lower triangle therefore ignored.
        continue;
    else
        currentIndex = currentIndex + 1;
    end
    
    jacobianM2rhoTmp(currentIndex, :) = PCx(k1,1:end) .* PCx(k2,1:end);
    
    [permutions1] = InnerCreatePermution(projection(k1,:));
    [permutions2] = InnerCreatePermution(projection(k2,:));
     
    PCD1 = transpose(projection(k2,:) * CxDrho);
    PCD2 = projection(k1,:) * CxDrho;

    
    part1 = permutions1 * PCD1;
    part2 = PCD2 * transpose(permutions2);
    jacobianM2XTmp(currentIndex,:) = part1 + transpose(part2);
    
end 

jacobianM2 = [jacobianM2rhoTmp, jacobianM2XTmp];

% Assamble
jacobian = (1 - pOutlier) * [jacobianM1; jacobianM2];

end


%% permutatin matrix function
% Create a matrix, which each row is shift of the vector
function [permutions] = InnerCreatePermution(rowVector)
 if(size(rowVector,1) ~= 1)
    rowVector = transpose(rowVector);
 end
 lengthVector = size(rowVector, 2);
 permutions = zeros(lengthVector, lengthVector);
 
 for index = 1 : lengthVector
     if ( index ~= 1)
        permutions(index,:) = [rowVector(1,index:end), rowVector(1,1:index - 1)];
    else
        permutions(index,:) = rowVector(1,:);
     end
 end 
end