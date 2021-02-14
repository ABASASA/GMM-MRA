%% add paths
clear;
clc;
AddPaths();
%% params
% Ground Truth input
L = 15; %Signal length

% Observations input
numberOfObservations = 100000; % Number of observations
sigmaArray = logspace(-2,0.5, 15); % sigma values to check
sigmaDiag = ones(L, 1); % homogenous noise case


% Comparassion input
numberOfStartingPoints = 2; % number of intial guesses
numberRepeats = 20; % This number represent how many signals & distributions (rhos) the code will compare).

% Optimization input
maxNumOfIterations = inf; % number of iteration in each GMM optimization
numberOfStepsGMM = 1; % number of steps in the GMM

% Paramter to each W
pOutlierStandard = 0; % outliers precent

k = 10;
projectionFunctionstandard =@() eye(L,L);%(k, L); % projection matrix

sigmaDiagStandard = ones(L, 1); % homogenous noise case
sigmaDiagHetro = 1:L; % homogenous noise case



% Define saving fig paramter
savingPath = 'Graphs/Compare_Geo_Dist/';

%% initialize data saveing objects
% Save SNR values
SNR = zeros(size(sigmaArray));
SNRHetro = zeros(size(sigmaArray));

geoStandated = zeros(numberRepeats, length(sigmaArray),2);
geoHetro = zeros(numberRepeats, length(sigmaArray),2);

%% Run on every noise level
for indexSigma = 1 : length(sigmaArray)
    % Compute sigma's covraice matrix
    sigmaStandard = (sigmaArray(indexSigma)^2) * diag(sigmaDiagStandard);
    sigmaHetro = (sigmaArray(indexSigma)^2) * diag(sigmaDiagHetro);

    currentSigmaOutlier = 10 * sigmaArray(indexSigma) *  eye(L,L);
    
    covOutlier = currentSigmaOutlier .^ 2; 
    % Compute SNR:
    SNR(indexSigma) = 1 /(sum(sigmaDiagStandard) * (sigmaArray(indexSigma)^2));
    SNRHetro(indexSigma) = 1 /(sum(sigmaDiagHetro) * (sigmaArray(indexSigma)^2));

    % create temperal data objection (for parfor) 
    
    geoStandatedTmp = zeros(numberRepeats, 2);
    geoHetroTmp = zeros(numberRepeats, 2);

    %% repets GMM
    parfor iRep = 1 : numberRepeats
        
        disp(num2str(iRep));
        %% Sample a ground truth (signal + shifts' distibution)
        
        [currentGroundTruth] = GenerateStartingPoints(L, 1);
        
        rho = currentGroundTruth(1:L);
        signal = currentGroundTruth(L + 1 : end);
        
        %% Create Projection 
        
        projectionStnadrd = projectionFunctionstandard();
        %%%%%%% Standard
        %% Create observations
        observations = GenerateObservations(signal, numberOfObservations,...
                            rho, sigmaStandard, projectionStnadrd, pOutlierStandard, currentSigmaOutlier);
                        
        %% Compute Estimations
        [M1Est, M2Est] = ComputeEmpricalMoments(observations);
        M2Est = ExtractUpperTriangleMatrixVectorize(M2Est);
        empricalMoment = [M1Est; M2Est];
        
        %% Define computation method to W
        momentFuction1by1 = @(theta, observations, sigma) ComputeMomentFucntion1By1(...
                  theta(1:L), theta((L+1) : end), observations, sigmaStandard,...
                  projectionStnadrd, pOutlierStandard, covOutlier);

        %% Compute W
        
        WStandart = ComputeW(currentGroundTruth , observations, sigmaStandard, momentFuction1by1);
        
        %%%%%%% Hetro

        %% Create observations
        observations = GenerateObservations(signal, numberOfObservations,...
                            rho, sigmaHetro, projectionStnadrd, pOutlierStandard, currentSigmaOutlier);
                        
        %% Compute Estimations
        [M1Est, M2Est] = ComputeEmpricalMoments(observations);
        M2Est = ExtractUpperTriangleMatrixVectorize(M2Est);
        empricalMoment = [M1Est; M2Est];
        
        %% Define computation method to W
        momentFuction1by1 = @(theta, observations, sigma) ComputeMomentFucntion1By1(...
                  theta(1:L), theta((L+1) : end), observations, sigmaHetro,...
                  projectionStnadrd, pOutlierStandard, covOutlier);

        %% Compute W
        
        WHetro = ComputeW(currentGroundTruth , observations, sigmaHetro, momentFuction1by1);
      
        
        %% Geo distance compute

        geoStandatedTmp(iRep,:) = [ComputeGeodeticDistanceFromI_uptoscalar(WStandart), cond(WStandart)];
        geoHetroTmp(iRep,:) = [ComputeGeodeticDistanceFromI_uptoscalar(WHetro), cond(WHetro)];

    end
    geoStandated(:, indexSigma,:) = geoStandatedTmp;
    geoHetro(:, indexSigma,:) = geoHetroTmp;


    disp(indexSigma);
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
end
%% Ask user if she/he wants to save the figures
disp(['Hey, You are running my code. Your Saveing Path is ', savingPath]);
disp('If you want to save the figure Press: 1 else 0');
toSave = input('Enter 1 or 0: ');
if toSave
    mkdir(savingPath)
    save([savingPath,'data.mat']);
end
%% Plot
fig = figure;
semilogx(SNR, mean(squeeze(geoStandated(:,:,1)),1), '*-', SNRHetro, mean(squeeze(geoHetro(:,:,1)),1), 's-');
legend({'Homoscedastic', 'Heteroscedastic'}, 'Location', 'northwest')
xlabel('SNR', 'fontsize', 12, 'fontweight','bold')
ylabel('Mean Geo. Distance', 'fontsize', 12, 'fontweight','bold')
xlim([0.01, 100])
fileName = 'Mean_Geo_Dist';
if toSave
    saveas(fig,[savingPath, fileName, '.fig']);
    saveas(fig,[savingPath, fileName, '.jpg']);
end

