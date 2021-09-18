function [jacobian] = ComputeJacobianM3MRA(rho, signal,projection, pOutlier, sigmaScalar )

[jacobianM1M2] = ComputeJacobianMRA(rho, signal, projection, pOutlier);

L =  length(signal);

[indeces] = ChooceIndecesM3(L);
numInd = size(indeces,1);

Csignal = circulant(signal); % circulant matrix of the signal

jacobianM3XTmp = zeros(numInd, L); 
jacobianM3RhoTmp = zeros(numInd, L); 


for iInd  = 1  : numInd
    
    k1  = indeces(iInd, 1);
    k2  = indeces(iInd, 2);
    k3  = indeces(iInd, 3);
    
    %% Rho
    xtmp1 = Csignal(k1, :);
    xtmp2 = Csignal(k2, :);
    xtmp3 = Csignal(k3, :);
    if (k1 == k2 & k2 == k3)
        bias = 3 * sigmaScalar^2 * xtmp1;
    elseif (k1 == k2)
        bias = sigmaScalar^2 * xtmp3;
    elseif (k2 == k3)
        bias = sigmaScalar^2 * xtmp1;
    elseif (k1 == k3)
        bias = sigmaScalar^2 * xtmp2;
    else
        bias = 0;
    end
    jacobianM3RhoTmp(iInd, :) =  (xtmp1 .* xtmp2 .* xtmp3 + bias);
    
    %% signal
    
    for iSignal = 1 : L
       
        j1 = mod(k1 - (iSignal),L) + 1;
        j2 = mod(k2 - (iSignal),L) + 1;
        j3 = mod(k3 - (iSignal),L) + 1;
        
        base = rho(j1) * signal(mod(k2 - (j1),L) + 1) * signal(mod(k3 - (j1),L) + 1) + ...
            rho(j2) * signal(mod(k1 - (j2),L) + 1) * signal(mod(k3 - (j2),L) + 1) +...
            rho(j3) * signal(mod(k1 - (j3),L) + 1) * signal(mod(k2 - (j3),L) + 1);
        
        if (k1 == k2 & k2 == k3)
            bias = 3 * sigmaScalar^2 * rho(j1);
        elseif (k1 == k2)
            bias = sigmaScalar^2 * rho(j3);
        elseif (k2 == k3)
            bias = sigmaScalar^2 * rho(j1);
        elseif (k1 == k3)
            bias = sigmaScalar^2 * rho(j2);
        else
            bias = 0;
        end
        
        jacobianM3XTmp(iInd, iSignal) = base + bias;
    end
    
end

jacobianM3 = [jacobianM3RhoTmp, jacobianM3XTmp];
jacobian = [jacobianM1M2;jacobianM3];

end