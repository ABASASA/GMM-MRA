function [indeces] = ChooceIndecesM3(L)

indeces = zeros(L^3, 3); % Full
counter = 1;
for k1 = 1 : L
    for k2 = k1 : L
        for k3 = k2 : L
            
            indeces(counter ,:) = [k1,k2,k3];
            counter = counter +1;
        end
    end
end
indeces = indeces(1: counter-1,:);
% indeces = [1,1,1];
end