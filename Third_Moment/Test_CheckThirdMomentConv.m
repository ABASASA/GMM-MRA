L = 7;
Ns = floor(logspace(2.5,5.2,8));

sigmaScalar = 0.5 ;
repNum = 10;

distances = inf(length(Ns), repNum);


sigmaDiag = ones(L, 1); % homogenous noise case
sigma = (sigmaScalar^2) * diag(sigmaDiag);

pOutlier = 0; % outliers precent
sigmaOutlier = 10 * eye(L,L);
projection = eye(L,L);


for  iN = 1 : length(Ns)
   
    N = Ns(iN);
    for iRep = 1 : repNum
        
        [currentGroundTruth] = GenerateStartingPoints(L, 1);
        
        rho = currentGroundTruth(1:L);
        signal = currentGroundTruth(L + 1 : end);
        
                
        %% Create observations
        observations = GenerateObservations(signal, N,...
                            rho, sigma, projection, pOutlier, sigmaOutlier);
        
        [indecesM3] = ChooceIndecesM3(L);
        [M3Ana] = ComuteM3Analytical(rho, signal, sigmaScalar,indecesM3);
        
        [M3Emp] = ComuteM3Empric(observations,indecesM3);
        
        
        distances(iN, iRep ) = norm(M3Ana(:) - M3Emp(:));
    end
    iN
end

res = mean(distances,2);
%% display
figure;
loglog(Ns, res, '*--')