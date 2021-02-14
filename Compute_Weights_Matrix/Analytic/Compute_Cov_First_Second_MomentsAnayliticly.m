%% Info
% Compute the covriance between each element of the first and second
% moment.


function [cov] = Compute_Cov_First_Second_MomentsAnayliticly(signal, rho, sigma, L)

%% Initlaize inces for sets
indeces = 1 : L;
A2 = [repmat(indeces, 1, L)]';
A1 = repmat(indeces, L, 1);
A1 = A1(:);

indeces = [A1, A2];

BL = size(indeces, 1);

%% Compute the moments
[M1, M2] = ComputeAnalyticMoments(rho, signal, sigma^2 * eye(L,L),...
                                eye(L,L), 0, eye(L,L));
estimateM2 = ExtractUpperTriangleMatrixVectorize(M2);

cov = zeros(length(estimateM2), L);

Cx  = circulant(signal);  % compute the circulent matrix of X
%% Compute
currentIndexR1 = 0;

for r1 = 1 : BL % Second Moment
    tmpset = indeces(r1,:);
    if (tmpset(1) < tmpset(2))
        continue;
    end
    currentIndexR1 = currentIndexR1 +1;
    for r2 = 1 : L % First Moment

        coeff = 0;
        set = [indeces(r1,:), r2];
        [C, ~, ic] = unique(set);
        %% First compute the forth moment 
        moment3 = rho' .* Cx(set(1), :) .* Cx(set(2), :) .* Cx(set(3), :);
        moment3 = sum(moment3);
        %% separte into cases
        if (length(C) == 3) % all the indeces are diffrent
            value = 0;
            % do nothing - only the 4-th moment
        elseif (length(C) == 2) % one pair
            value = CaseOfOnePair(C, ic, sigma,  Cx, rho');
        elseif (length(C) == 1) % Quadruple
            value = CaseOfTriplate(C, ic, sigma,  Cx, rho');
        end
        
        coeff = value + moment3;

         cov(currentIndexR1,r2) = coeff;
    end
end

%%

cov = cov - estimateM2 * M1';
end


function [value] = CaseOfOnePair(C,ic, sigma, Cx, rho)
    %% first need to find the pair
    indexes = 1:2;
    for index = 1 : 2
        if sum(ic == index) == 2 % Found the pair
            relevantIndexesInC = C(indexes(indexes ~= index)); % Found the other 2 indexes
        end
    end
    value = (sigma^2) * sum(rho .* Cx(relevantIndexesInC(1),:));
end

function [value] = CaseOfTriplate(C,ic, sigma, Cx, rho) 
    value = 3 * (sigma^2) *sum(rho .* Cx(C(1),:));
end

