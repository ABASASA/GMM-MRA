function [M1, M2] = ComputeMomentsEsimateL1GeometriclMedian(observations, sigma, projection, batchSize)

N = size(observations, 2);
L = size(observations, 1);

numbrtOfBatches = ceil(N / batchSize);
EachBatchResultsM1 = zeros(L, numbrtOfBatches);
EachBatchResultsM2 = zeros(L^2, numbrtOfBatches);

for iBatch = 1 : numbrtOfBatches
    if  iBatch == numbrtOfBatches
        currentObservations = observations(:,...
                            (iBatch-1) * batchSize + 1: end);
    else
        currentObservations = observations(:,...
                            (iBatch-1) * batchSize + 1: iBatch * batchSize);
    end
    Ncurrent = size(currentObservations, 2);
    dataM1 = struct;
    dataM1.Data = currentObservations';
    dataM1.RelTol = 0.0001;
%     EachBatchResultsM1(:,iBatch) = geometric_median(currentObservations);
    tmpDataM1 = Weiszfeld(dataM1);
    EachBatchResultsM1(:,iBatch) = tmpDataM1.xMedian;
    M2s = zeros(L^2, Ncurrent);
    M2 = zeros(L, L);
    %% Compute for each sample 
    for i = 1 : Ncurrent
        tmpVar = observations(:,i) * transpose(currentObservations(:, i));
        M2s(:,i) =tmpVar(:);
    end
    %% Add projection and un-bias
%     EachBatchResultsM2(:,iBatch) = geometric_median(M2s);
    dataM2 = struct;
    dataM2.Data = M2s';
    dataM2.RelTol = 0.0001;
    tmpDataM2 = Weiszfeld(dataM2);
    EachBatchResultsM2(:,iBatch) = tmpDataM2.xMedian;

end
if (numbrtOfBatches > 1)
    dataM1.Data = EachBatchResultsM1';
    dataM2.Data = EachBatchResultsM2';
    dataM1.RelTol = 0.0001;
    dataM2.RelTol = 0.0001;

    % M1 = geometric_median(EachBatchResultsM1);
    % M2 = geometric_median(EachBatchResultsM2);
    tmpDataM1 = Weiszfeld(dataM1);
    tmpDataM2 = Weiszfeld(dataM2);
    M1 = tmpDataM1.xMedian;
    M2 = tmpDataM2.xMedian;
else
    M1 = EachBatchResultsM1;
    M2 = EachBatchResultsM2;
end
    
M2 = reshape(M2,[L,L])  -  projection * sigma *  transpose(projection) ; 

end


function output_structure = Weiszfeld(input_structure)
%
% This function numerically calculates the geometric mean of a
% N-Dimensional set of points using the Wieszfeld's algorithm
%
% INPUT: Structure
% -- [REQUIRED] input_strucuture.Data = Input data matrix. Each row is a point, each
% column a dimension
% -- input_structure.RelTol = Relative tolerance for stopping the search.
% Default: 0.001
% -- input_structure.x0 = A vector with the initial point. If not provided,
% it is automatically calcualted based on the centroid of the original
% series of points
%% Default values
RelTolDefault = 0.001 ;
expectedIterations = 20 ;
%% Reading inputs
% Reading the data matrix 
data = input_structure.Data ;
% Check if the RelTol field is provided. Otherwise, use the default value
if any(strcmp(fields(input_structure),'RelTol'))
    relTol = input_structure.RelTol ;
else
    relTol = RelTolDefault ;
end
% Check if a starting point is provided. Otherwise, calculate it
if any(strcmp(fields(input_structure),'x0'))
    x0 = input_structure.x0 ;
else
    x0 = mean(data,1) ;
end
%% Calculating some useful parameters
[nPoints, nDimensions] = size(data) ;
% Initialize the relative difference
eps = 1 ;
counter = 0 ;
% Initialize the matrix storing all iterations. 
xTemp = NaN([expectedIterations , nDimensions]) ;
xTemp(1,:) = x0 ;
%% Iterations
while eps > relTol
    counter = counter + 1 ;
    % Chaged here from expilciy to inexpilciy (bsxfun)
    weights = sum((bsxfun(@minus,data,xTemp(counter,:))).^2, 2).^-0.5 ;
    repMatWeights = repmat(weights,[1,nDimensions]);
    temp = sum(data .* (repMatWeights .* ones(nPoints, nDimensions)), 1) / sum(weights) ;
    xTemp(counter+1, :) = temp ;
    eps = (sum((xTemp(counter+1,:) - xTemp(counter,:)).^2))^0.5 ;
end
%% Post compute
% Compute the difference at the last computation
err = sum(sum((bsxfun(@minus,data,xTemp(counter,:))).^2,2).^0.5) / nPoints;
output_structure.xMedian = xTemp(counter+1,:) ;
output_structure.err = err ;
output_strucutre.tol = eps ;
end