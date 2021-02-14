%% Info
% This function compute the variance of the seecond moments analyticly.
% Asaf Abas 4.9.20
function [W] = ComputeVarSecondMomentsAnayliticly(signal, rho, sigma, L)

%% Initlaize inces for sets
indeces = 1 : L;
A2 = [repmat(indeces, 1, L)]';
A1 = repmat(indeces, L, 1);
A1 = A1(:);

indeces = [A1, A2];

BL = size(indeces, 1);

%% Compute the moments
[~, M2] = ComputeAnalyticMoments(rho, signal, sigma.^2 * eye(L,L),...
                                eye(L,L), 0, eye(L,L));
estimateM2 = ExtractUpperTriangleMatrixVectorize(M2);

W = zeros(length(estimateM2), length(estimateM2));

Cx  = circulant(signal);  % compute the circulent matrix of X
%% Compute  
currentIndexR1 = 0;
for r1 = 1 : BL
    tmpset = indeces(r1,:);
    if (tmpset(1) < tmpset(2))
        continue;
    end
    currentIndexR1 = currentIndexR1 +1;
    currentIndexR2 = 0;
    for r2 = 1 : BL

        coeff = 0;
        set = [indeces(r1,:), indeces(r2,:)];
        if (set(3) < set(4))
            continue;
        end
        currentIndexR2 = currentIndexR2 + 1;
        [C, ~, ic] = unique(set);
        %% First compute the forth moment 
        moment4 = rho' .* Cx(set(1), :) .* Cx(set(2), :) .* Cx(set(3), :) .* ...
                         Cx(set(4), :);
        moment4 = sum(moment4);
        %% separte into cases
        if (length(C) == 4) % all the indeces are diffrent
            value = 0;
            % do notthing - only the 4-th moment
        elseif (length(C) == 3) % one pair
            value = CaseOfOnePair(C, ic, sigma, Cx, rho');
            
        elseif (length(C) == 2) % 2 Pairs or Triplate
            if (sum(ic == 1) == 2) % This is a 2 pairs
                value = CaseOfTwoPair(C, ic, sigma, Cx, rho');
            else % This case is a triplare
                value = CaseOfTriplate(C, ic, sigma, Cx, rho');
            end
        elseif (length(C) == 1) % Quadruple
            value = CaseOfQuadruple(C, ic, sigma, Cx, rho');
        end
        
        coeff = value + moment4;

         W(currentIndexR1, currentIndexR2) = coeff;
    end
end

W = W -  estimateM2 * estimateM2';
end


function [value] = CaseOfOnePair(C,ic, sigma, Cx, rho)
    %% first need to find the pair
    indexes = 1:3;
    for index = 1 : 3
        if sum(ic == index) == 2 % Found the pair
            relevantIndexesInC = C(indexes(indexes ~= index)); % Found the other 2 indexes
        end
    end
    value = (sigma^2) * sum(rho .* Cx(relevantIndexesInC(1),:) .*...
                            Cx(relevantIndexesInC(2),:));
end

function [value] = CaseOfTwoPair(C,ic, sigma, Cx, rho)
    firstPair = (sigma^2) * sum(rho .* Cx(C(1),:) .*...
                                Cx(C(1),:) );
    secondPair = (sigma^2) * sum(rho .* Cx(C(2),:) .*...
                                 Cx(C(2),:) );
    value = firstPair + secondPair + (sigma^4);
end


function [value] = CaseOfTriplate(C,ic, sigma, Cx, rho) 

    value = 3 * (sigma^2) * sum(rho .* Cx(C(1),:) .*...
                                 Cx(C(2),:) );
end

function [value] = CaseOfQuadruple(C,ic, sigma, Cx, rho)
    partOne = 6 * (sigma^2) * sum(rho .* Cx(C(1),:) .*...
                                 Cx(C(1),:));
    partTwo = 3 * (sigma^4);
    value = partOne + partTwo;
end
