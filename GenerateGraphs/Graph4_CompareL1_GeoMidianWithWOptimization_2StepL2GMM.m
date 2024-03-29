%% add paths
clear;
clc;
AddPaths();
%% params
% Ground Truth input
L = 15; %Signal length

% Observations input
numberOfObservations = 100000; % Number of observations
sigmaArray = logspace(-1.8,0.5,8);%logspace(-2,1.2,10); % sigma values to check
sigmaDiag = ones(L,1); % homogenous noise case
pOutlier = 0.2; % outliers precent

projectionFunction =@() eye(L,L);%(k, L); % projection matrix

% Comparassion input
numberOfStartingPoints = 2; % number of intial guesses
numberRepeats = 80; % This number represent how many signals & distributions (rhos) the code will compare).

% Optimization input
maxNumOfIterations = inf; % number of iteration in each GMM optimization
numberOfStepsGMM = 2; % number of steps in the GMM

% Define saving fig paramter
savingPath = 'Graphs/L1_GeoMid_WithWOptimization_vs_2-Step_L2_GMM_outliers_02-countFlags/';

%% initialize data saveing objects
% Save SNR values
SNR = zeros(size(sigmaArray));

relativeErrorLS = zeros(length(sigmaArray), numberRepeats, 4); % one row for signal and one for distibution.
relativeErrorGMM = zeros(length(sigmaArray), numberRepeats, 4);
flags = inf([length(sigmaArray),numberRepeats, numberOfStartingPoints, 2]);

%% Run on every noise level
for indexSigma = 1 : length(sigmaArray)
    % Compute sigma's covraice matrix
    sigma = (sigmaArray(indexSigma)^2) * diag(sigmaDiag);
    sigmaOutliers = 10 * sigmaArray(indexSigma) * eye(L,L);
    covOutliers = (sigmaOutliers).^2 * eye(L,L);
    % Compute SNR:
    SNR(indexSigma) = 1 /(sum(sigmaDiag) * (sigmaArray(indexSigma)^2));
    
    % create temperal data objection (for parfor) 
    
    tmprelativeErrorLS = zeros(numberRepeats, 4);
    tmprelativeErrorGMM = zeros(numberRepeats, 4);
    tmpFlagsLS = inf([numberRepeats, numberOfStartingPoints]);
    tmpFlagsGMM = inf([numberRepeats, numberOfStartingPoints]);

    %% repets GMM
    parfor iRep = 1 : numberRepeats
        
        disp(num2str(iRep));
        %% Sample a ground truth (signal + shifts' distibution)
        
        [currentGroundTruth] = GenerateStartingPoints(L, 1);
        
        rho = currentGroundTruth(1:L);
        signal = currentGroundTruth(L + 1 : end);
        
        %% Create Projection 
        
        projection = projectionFunction();
        
        %% Create observations
        observations = GenerateObservations(signal, numberOfObservations,...
                            rho, sigma, projection, pOutlier, sigmaOutliers);
                        
        %% Compute Estimations
        [M1Est, M2Est] = ComputeEmpricalMoments(observations);
        M2Est = ExtractUpperTriangleMatrixVectorize(M2Est);
        empricalMomentMean = [M1Est; M2Est];
        
        %% Compute Estimation GeometriclMedian
        [M1Est, M2Est] = ComputeMomentsEsimateL1GeometriclMedian(observations, sigma, projection, numberOfObservations + 1);
        M2Est = ExtractUpperTriangleMatrixVectorize(M2Est);
        empricalMomentGeoMedian = [M1Est; M2Est];
        %% Random a starting point for the optimization
        [startingPoints] = GenerateStartingPoints(L, numberOfStartingPoints);

        %% Define computation method to W
        momentFuction1by1 = @(theta, observations, sigma) ComputeMomentFucntion1By1(...
                  theta(1:L), theta((L+1) : end), observations, sigma,...
                  projection, pOutlier, covOutliers);
        momentFuction1by1Direct = @(theta, observations, sigma) ComputeMomentFucntion1By1Direct(...
                  theta(1:L), theta((L+1) : end), observations, sigma,...
                  projection, pOutlier, covOutliers);     
        %% Compute W
        
        WL1 = ones(length(empricalMomentMean),1);
        WL1 = WL1 ./ norm(WL1);
        WGMM = eye(length(empricalMomentMean), length(empricalMomentMean));
        WGMM = ComputeW(currentGroundTruth , observations, sigma, momentFuction1by1Direct);

        %% LS
        [estSignalLS, estRhoLS, infoLS, ~] = Iterative_GMM_L1_Optimization(startingPoints,...
                WL1, empricalMomentGeoMedian, sigma, projection, pOutlier,...
                covOutliers, maxNumOfIterations,...
                numberOfStepsGMM, 1);
        
        %% Compute relative Errors - LS
        relativeErrorSignalLS = RelativeErrorUpToShift(signal, estSignalLS);
        relativeErrorRhoLS = RelativeErrorUpToShift(rho, estRhoLS);
        numberOfIterLS = infoLS.iterations;
        cpuTimeLS = infoLS.CPUTime;
        tmprelativeErrorLS(iRep,:) =  [relativeErrorSignalLS,...
            relativeErrorRhoLS, numberOfIterLS, cpuTimeLS];
         disp([num2str(iRep) ': end L1']);

        %% Opt GMM
        [estSignalGMM, estRhoGMM, infoGMM, ~] = ComputeIterativeGMMviaMatlab(...
                startingPoints, WGMM, observations, empricalMomentMean, sigma,...
                momentFuction1by1, projection, pOutlier, covOutliers,...
                maxNumOfIterations, 1, 1);
        %% Compute relative Errors - GMM
        relativeErrorSignalGMM = RelativeErrorUpToShift(signal, estSignalGMM);
        relativeErrorRhoGMM = RelativeErrorUpToShift(rho, estRhoGMM);
        numberOfIterGMM = infoGMM.iterations;
        cpuTimeGMM = infoGMM.CPUTime;%singularValueDistanceToEye(WLS);
        tmprelativeErrorGMM(iRep,:) =  [relativeErrorSignalGMM,...
            relativeErrorRhoGMM, numberOfIterGMM, cpuTimeGMM];
        disp([num2str(iRep) ': end GMM']);
        tmpFlagsLS(iRep,:)= infoLS.Flags ;
        tmpFlagsGMM(iRep,:)= infoGMM.Flags ;

       % tmpFlags( (iRep-1)* numberRepeats : iRep* numberRepeats, 1)= infoGMM.Flags;
        %tmpFlags( (iRep-1)* numberRepeats : iRep* numberRepeats, 2)= infoGMM.Flags;
    end
    flags(indexSigma,:,:,1)      = tmpFlagsLS;
    flags(indexSigma,:,:,2)      = tmpFlagsGMM;

    relativeErrorLS(indexSigma,:,:) = tmprelativeErrorLS;
    relativeErrorGMM(indexSigma,:,:) = tmprelativeErrorGMM;
    
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp(['Yo man, im in index....' num2str(indexSigma) '  out of...' num2str(length(sigmaArray))]);
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
end

%% Compare - Relative Error (1 - signal, 2 - rho)
% Ask user if she/he wants to save the figures
disp(['Hey, You are running my code. Your Saveing Path is ', savingPath]);
disp('If you want to save the figure Press: 1 else 0');
% toSave = input('Enter 1 or 0: ');
toSave = 1;
if(toSave)
    mkdir(savingPath);
    save([savingPath,'data.mat']);
end
for ii = 1 : 4
    meanrelativeErrorLS = zeros(length(sigmaArray),1);
    meanrelativeErrorOpt = zeros(length(sigmaArray),1);

    
    stdrelativeErrorLS = zeros(length(sigmaArray),1);
    stdrelativeErrorOpt = zeros(length(sigmaArray),1);
    proportions = zeros(length(sigmaArray),numberRepeats);
    %% Compute mean & variance for each sigma
    for indexSigma = 1 : length(sigmaArray)
%         meanrelativeErrorLS(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
%         meanrelativeErrorOpt(indexSigma) = mean(squeeze(relativeErrorGMM(indexSigma,:,ii)));
%         
%         stdrelativeErrorLS(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
%         stdrelativeErrorOpt(indexSigma) = std(squeeze(relativeErrorGMM(indexSigma,:,ii)));
%         % Ratio
%         proportions(indexSigma,:) = (squeeze(relativeErrorLS(indexSigma,:,ii))./...
%                                     squeeze(relativeErrorGMM(indexSigma,:,ii)));
                                
        meanrelativeErrorLS(indexSigma) = mean(squeeze(relativeErrorGMM(indexSigma,:,ii)));
        meanrelativeErrorOpt(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
        
        stdrelativeErrorLS(indexSigma) = std(squeeze(relativeErrorGMM(indexSigma,:,ii)));
        stdrelativeErrorOpt(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
        % Ratio
        proportions(indexSigma,:) = (squeeze(relativeErrorGMM(indexSigma,:,ii))./...
                                    squeeze(relativeErrorLS(indexSigma,:,ii)));
    end
    

    %% Plot histogram - Relative Error
    fig = figure();
    %% Mean
    subplot(2,1,1);
    semilogx(SNR, meanrelativeErrorLS, 'r*-',...
        SNR, meanrelativeErrorOpt, 'b*-');

    xlabel('SNR');
    legend( 'GMM', 'L1-Geo. Median', 'Location','southwest');
    if (ii == 1)
         title(['Signal Mean - Relative Error']);
    elseif (ii == 2)
        title(['\rho Mean - Relative Error']);
    elseif (ii == 3)
        title(['Mean Number Of Iterations']);
    elseif (ii == 4)
        title(['Mean CPU Time']);
    end
    %% Variance
    subplot(2,1,2);
    semilogx(SNR, stdrelativeErrorLS, 'r*-', ...
        SNR, stdrelativeErrorOpt, 'b*-');
    
    xlabel('SNR');
    legend( 'GMM', 'L1-Geo. Median', 'Location','southwest');
    if (ii == 1)
         title(['Signal std - Relative Error']);
         if (toSave)
            fileName = 'Signal  - Relative Error';
            mkdir(savingPath);
            saveas(fig,[savingPath, fileName, '.fig']);
            saveas(fig,[savingPath, fileName, '.jpg']);
        end
    elseif (ii == 2)
        title(['\rho std - Relative Error']);
        if (toSave)
            fileName = 'Rho - Relative Error';
            saveas(fig,[savingPath, fileName, '.fig']);
            saveas(fig,[savingPath, fileName, '.jpg']);
        end
    elseif (ii == 3)
        title(['std Number Of Iterations']);
        
        if (toSave)
            fileName = 'Number Of Iterations';
            saveas(fig,[savingPath, fileName, '.fig']);
            saveas(fig,[savingPath, fileName, '.jpg']);
        end
    elseif (ii == 4)
        title(['std CPU Time']);
        
        if (toSave)
            fileName = 'CPU Time';
            saveas(fig,[savingPath, fileName, '.fig']);
            saveas(fig,[savingPath, fileName, '.jpg']);
        end
    end
    
     fig = figure;
   % [fig] = BoxPlotAsaf(fig, SNR(end:-1:1), proportions(end:-1:1,end:-1:1), 'b*-');
    [fig] = BoxPlotAsaf(fig, SNR, proportions, 'b*-');
    [indexMeanBigger1] = find(meanrelativeErrorOpt >= 1,1);
    hold on;
    plot(1: size(proportions,1), ones(size(proportions,1),1),'k--');
    hold on;
    if ~isempty(indexMeanBigger1)
        plot((indexMeanBigger1) * ones(2,1), [max(proportions(:)); min(proportions(:))],'g');
    end
    xlabel('\sigma');
    xlabel('SNR');
    ylabel('Ratio');
    if (ii == 1)
        title('Ratio Signal Relative Error - L2-GMM / L1-Geo. Median');
    elseif (ii == 2)
        title('Ratio \rho Relative Error  -  L2-GMM / L1-Geo. Median');
    elseif (ii == 3)
        title('Ratio #iter - L2-GMM / L1-Geo. Median');
    elseif (ii == 4)
        title('Ratio CPU Time -  L2-GMM / L1-Geo. Median');
    end
    if (toSave)
        if (ii == 1)
            fileName = 'Ratio_Relative_Error';
        elseif (ii == 2)
            fileName = 'Ratio_Disribution_Relative_Error';
        elseif (ii == 3)
            fileName = 'Ratio_Number_Iterations';
        elseif (ii == 4)
            fileName = 'Ratio_CPU_Time';
        end
        saveas(fig,[savingPath, fileName, '.fig']);
        saveas(fig,[savingPath, fileName, '.jpg']);
    end

end
MakeGraphPretty_Outliers;