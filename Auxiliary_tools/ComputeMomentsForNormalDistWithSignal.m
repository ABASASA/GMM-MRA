function [moments] = ComputeMomentsForNormalDistWithSignal(sigma,x)

moments = zeros(length(x),6);
for  i = 1 : length(x)
       moments(i,:) = First6MomentsOfNormalDist(x(i) - mean(x), sigma); 
end


end
