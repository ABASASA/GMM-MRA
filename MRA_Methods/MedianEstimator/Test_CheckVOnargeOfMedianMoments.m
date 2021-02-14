%% User Inputs
Ns = floor(logspace(3,6,9));
% Ns = Ns(3)
repets = 150;
L = 13;
sigmaNoise = 0.1;
pOutlier = 0.25; % p is 100*p% to get an outlier
pEstiamteOutlier = 0.2;
projection = eye(L, L);
GMBatchSize = Ns(end)+1;
%% Saveing data
savingPath = 'Graphs/Test_MRA_Outliers_FixedAnalyticEstimator030620/Sigma_01_pOutTrue_025_PEst_02/';
%% Init
sigma = sigmaNoise * eye(L,L);
sigmaOutlier = 10 * sigmaNoise;
CovOutlier = sigmaOutlier^2 * eye(L,L);
M1ErrorMarg = zeros(length(Ns), repets);
M2ErrorMarg = zeros(length(Ns), repets);

M1ErrorGeo = zeros(length(Ns), repets);
M2ErrorGeo = zeros(length(Ns), repets);

M1ErrorMean = zeros(length(Ns), repets);
M2ErrorMean = zeros(length(Ns), repets);

for iN = 1 : length(Ns)
    N = Ns(iN);
    M1ErrorTmpMarg = zeros(repets, 1);
    M2ErrorTmpMarg = zeros(repets, 1);

    M1ErrorTmpGeo = zeros(repets, 1);
    M2ErrorTmpGeo = zeros(repets, 1);

    M1ErrorTmpMean = zeros(repets, 1);
    M2ErrorTmpMean = zeros(repets, 1);
    %% Loop
    parfor iRep = 1 : repets

        %% Random a signal
        signal = randn(L,1);
        signal = signal / norm(signal);
        %% Create Rho
        rho = rand(L,1); % posative number between [0,1]
        rho = rho ./ sum(rho); % posative number between [0,1] and sum is 1
        %% Create observations
        observations = RandSignalbyRhoWithOutliers(signal, N, rho,...
                                            sigma, projection, pOutlier,sigmaOutlier);
        %% Moments From estimations
        % Marginal medain
        [M1EstMarg, M2EstMarg] = ComputeMomentsEsimateL1MatginalMedian(observations, sigma, projection, N);
        M2EstMarg = ExtractUpperTriangleMatrixVectorize(M2EstMarg);
        
        % Geometric Meedian
        [M1EstGeo, M2EstGeo] = ComputeMomentsEsimateL1GeometriclMedian(observations, sigma, projection, N, GMBatchSize);
        M2EstGeo = ExtractUpperTriangleMatrixVectorize(M2EstGeo);
%         estimations = [M1Est; M2Est];
        % Mean estimator
        [M1Mean, M2Mean] = ComputeMomentsEsimateOutliers(observations, sigma,...
                                        projection, N, pEstiamteOutlier);
        M2Mean = ExtractUpperTriangleMatrixVectorize(M2Mean);
        %% True mometns

        [M1, M2] = ComputeMomentsAccurateOutliers(rho, signal, projection,...
                                                pEstiamteOutlier,CovOutlier);
        M2 = ExtractUpperTriangleMatrixVectorize(M2);

        %% Save data
        M1ErrorTmpMarg(iRep) = norm(M1EstMarg - M1,2) / length(M1);
        M2ErrorTmpMarg(iRep) = norm(M2EstMarg - M2,2) / length(M2);
        
        M1ErrorTmpGeo(iRep) = norm(M1EstGeo - M1,2) / length(M1);
        M2ErrorTmpGeo(iRep) = norm(M2EstGeo - M2,2) / length(M2);
        
        M1ErrorTmpMean(iRep) = norm(M1Mean - M1,2) / length(M1);
        M2ErrorTmpMean(iRep) = norm(M2Mean - M2,2) / length(M2);
        
    end
    M1ErrorMarg(iN,:) = M1ErrorTmpMarg;
    M2ErrorMarg(iN,:) = M2ErrorTmpMarg;
    
    M1ErrorGeo(iN,:) = M1ErrorTmpGeo;
    M2ErrorGeo(iN,:) = M2ErrorTmpGeo;
    
    M1ErrorMean(iN,:) = M1ErrorTmpMean;
    M2ErrorMean(iN,:) = M2ErrorTmpMean;
    disp(['Finished ' num2str(iN) ' out of ' num2str(length(Ns))]);
end
%%
meanM1ErrorMarg = mean(M1ErrorMarg,2);
meanM2ErrorMarg = mean(M2ErrorMarg,2);

meanM1ErrorGeo = mean(M1ErrorGeo,2);
meanM2ErrorGeo = mean(M2ErrorGeo,2);

meanM1ErrorMean = mean(M1ErrorMean,2);
meanM2ErrorMean = mean(M2ErrorMean,2);


medianM1ErrorMarg = median(M1ErrorMarg,2);
medianM2ErrorMarg = median(M2ErrorMarg,2);

medianM1ErrorGeo = median(M1ErrorGeo,2);
medianM2ErrorGeo = median(M2ErrorGeo,2);

medianM1ErrorMean = median(M1ErrorMean,2);
medianM2ErrorMean = median(M2ErrorMean,2);

%% plot
disp(['Hey, You are running my code. Your Saveing Path is ', savingPath]);
disp('If you want to save the figure Press: 1 else 0');
toSave = input('Enter 1 or 0: ');
if(toSave)
    mkdir(savingPath);
    save([savingPath,'data.mat']);
end

% Marginal Median

fig1 = figure;
subplot(2,1,1);
loglog(Ns,meanM1ErrorMarg ,'b*--', Ns, medianM1ErrorMarg,'r*--');
title('Marginal - Error in M1');
xlabel('#Observations')

subplot(2,1,2);
loglog(Ns,meanM2ErrorMarg ,'b*--', Ns, medianM2ErrorMarg,'r*--');
title('Marginal - Error in M2');
legend('Mean Error','Median Error','Location','southwest');
xlabel('#Observations')

% Geo Median
fig2 = figure;
subplot(2,1,1);
loglog(Ns,meanM1ErrorGeo ,'b*--', Ns, medianM1ErrorGeo,'r*--');
title('Geo - Error in M1');
xlabel('#Observations')

fig3 = subplot(2,1,2);
loglog(Ns,meanM2ErrorGeo, 'b*--', Ns, medianM2ErrorGeo,'r*--');
title('Geo - Error in M2');
legend('Mean Error','Median Error','Location','southwest');
xlabel('#Observations')

% Mean
figure;
subplot(2,1,1);
loglog(Ns,meanM1ErrorMean ,'b*--', Ns, medianM1ErrorMean,'r*--');
title('Mean - Error in M1');
xlabel('#Observations')

subplot(2,1,2);
loglog(Ns,meanM2ErrorMean, 'b*--', Ns, medianM2ErrorMean,'r*--');
title('Mean - Error in M2');
legend('Mean Error','Median Error','Location','southwest');
xlabel('#Observations')
if(toSave)
    fileName1 = 'MarginalMedianEstimator';
    saveas(fig1,[savingPath, fileName1, '.fig']);
    saveas(fig1,[savingPath, fileName1, '.jpg']);

    fileName2 = 'GeometricMedianEstimator';
    saveas(fig2,[savingPath, fileName2, '.fig']);
    saveas(fig2,[savingPath, fileName2, '.jpg']);

    fileName3 = 'MeanEstimator';
    saveas(fig3,[savingPath, fileName3, '.fig']);
    saveas(fig3,[savingPath, fileName3, '.jpg']);

end