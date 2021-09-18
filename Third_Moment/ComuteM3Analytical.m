function [M3] = ComuteM3Analytical(rho, signal, sigma, indeces)

L = length(signal); % Length
Csignal = circulant(signal); % circulant matrix of the signal

numInd = size(indeces,1);
M3 = zeros(numInd , 1);

for ind = 1 : numInd
    k1  = indeces(ind, 1);
    k2  = indeces(ind, 2);
    k3  = indeces(ind, 3);

    xtmp1 = Csignal(k1, :);
    xtmp2 = Csignal(k2, :);
    xtmp3 = Csignal(k3, :);
    if (k1 == k2 & k2 == k3)
        bias = 3 * sigma^2 * xtmp1;
    elseif (k1 == k2)
        bias = sigma^2 * xtmp3;
    elseif (k2 == k3)
        bias = sigma^2 * xtmp1;
    elseif (k1 == k3)
        bias = sigma^2 * xtmp2;
    else
        bias = 0;
    end
    M3(ind) = sum(transpose(rho) .* (xtmp1 .* xtmp2 .* xtmp3 + bias) );
end

end