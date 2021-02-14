%% add paths
clear;
clc;
AddPaths();
%% params
% Ground Truth input
L = 15; %Signal length

% Observations input
numberOfObservations = 100000; % Number of observations
sigmaArray = logspace(-2,1,8); % sigma values to check
sigmaDiag = ones(L, 1); % homogenous noise case
pOutlier = 0; % outliers precent
sigmaOutlier = 10 * eye(L,L);

projectionFunction =@() eye(L,L);%(k, L); % projection matrix

% Comparassion input
numberOfStartingPoints = 2; % number of intial guesses
numberRepeats = 100; % This number represent how many signals & distributions (rhos) the code will compare).

% Optimization input
maxNumOfIterations = inf; % number of iteration in each GMM optimization
numberOfStepsGMM = 1; % number of steps in the GMM
numberOfStepsGMM2Step = 2;
% Define saving fig paramter
savingPath = 'Graphs/Compare_LS_OptimalGMM_2StepGMM-StandardCase/';

%% initialize data saveing objects
% Save SNR values
SNR = zeros(size(sigmaArray));

relativeErrorLS = zeros(length(sigmaArray), numberRepeats, 4); % one row for signal and one for distibution.
relativeErrorGMMOptimal = zeros(length(sigmaArray), numberRepeats, 4);
relativeErrorGMM2Step = zeros(length(sigmaArray), numberRepeats, 4);


%% Run on every noise level
for indexSigma = 1 : length(sigmaArray)
    % Compute sigma's covraice matrix
    sigma = (sigmaArray(indexSigma)^2) * diag(sigmaDiag);
    sigmaOurlierCurrent =  sigmaArray(indexSigma) * sigmaOutlier;
    covOurlier = sigmaOurlierCurrent .^ 2;
    % Compute SNR:
    SNR(indexSigma) = 1 /(sum(sigmaDiag) * (sigmaArray(indexSigma)^2));
    
    % create temperal data objection (for parfor) 
    
    tmprelativeErrorLS = zeros(numberRepeats, 4);
    tmprelativeErrorGMMOptimal = zeros(numberRepeats, 4);
    tmprelativeErrorGMM2Step = zeros(numberRepeats, 4);
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
                            rho, sigma, projection, pOutlier, sigmaOurlierCurrent);
                        
        %% Compute Estimations
        [M1Est, M2Est] = ComputeEmpricalMoments(observations);
        M2Est = ExtractUpperTriangleMatrixVectorize(M2Est);
        empricalMoment = [M1Est; M2Est];
        
        %% Random a starting point for the optimization
        [startingPoints] = GenerateStartingPoints(L, numberOfStartingPoints);

        %% Define computation method to W
        momentFuction1by1 = @(theta, observations, sigma) ComputeMomentFucntion1By1(...
                  theta(1:L), theta((L+1) : end), observations, sigma,...
                  projection, pOutlier, covOurlier);

        %% Compute W
        
        [WGMM] = ComputeW(currentGroundTruth , observations, sigma, momentFuction1by1);
        WLS = eye(size(WGMM));

        %% LS
        [estSignalLS, estRhoLS, infoLS, ~] = ComputeIterativeGMMviaMatlab(startingPoints,...
                WLS, observations, empricalMoment, sigma,...
                momentFuction1by1, projection, pOutlier, covOurlier,...
                maxNumOfIterations, 1, 1);
        
        %% Compute relative Errors - LS
        relativeErrorSignalLS = RelativeErrorUpToShift(signal, estSignalLS);
        relativeErrorRhoLS = RelativeErrorUpToShift(rho, estRhoLS);
        numberOfIterLS = infoLS.iterations;
        cpuTimeLS = infoLS.CPUTime;
        tmprelativeErrorLS(iRep,:) =  [relativeErrorSignalLS,...
            relativeErrorRhoLS, numberOfIterLS, cpuTimeLS];
        %% Opt GMM
        [estSignalGMM, estRhoGMM, infoGMM, ~] = ComputeIterativeGMMviaMatlab(...
                startingPoints, WGMM, observations, empricalMoment, sigma,...
                momentFuction1by1, projection, pOutlier, covOurlier,...
                maxNumOfIterations, numberOfStepsGMM, 1);
        % Compute relative Errors - GMM
        relativeErrorSignalGMM = RelativeErrorUpToShift(signal, estSignalGMM);
        relativeErrorRhoGMM = RelativeErrorUpToShift(rho, estRhoGMM);
        numberOfIterGMM = infoGMM.iterations;
        cpuTimeGMM = infoGMM.CPUTime;%singularValueDistanceToEye(WLS);
        tmprelativeErrorGMMOptimal(iRep,:) =  [relativeErrorSignalGMM,...
            relativeErrorRhoGMM, numberOfIterGMM, cpuTimeGMM];
      %% 2-steps GMM
        [estSignalGMM, estRhoGMM, infoGMM, ~] = ComputeIterativeGMMviaMatlab(...
                startingPoints, WLS, observations, empricalMoment, sigma,...
                momentFuction1by1, projection, pOutlier, covOurlier,...
                maxNumOfIterations, numberOfStepsGMM2Step, 1);
        % Compute relative Errors - GMM
        relativeErrorSignalGMM = RelativeErrorUpToShift(signal, estSignalGMM);
        relativeErrorRhoGMM = RelativeErrorUpToShift(rho, estRhoGMM);
        numberOfIterGMM = infoGMM.iterations;
        cpuTimeGMM = infoGMM.CPUTime;%singularValueDistanceToEye(WLS);
        tmprelativeErrorGMM2Step(iRep,:) =  [relativeErrorSignalGMM,...
            relativeErrorRhoGMM, numberOfIterGMM, cpuTimeGMM];
    end

    relativeErrorLS(indexSigma,:,:) = tmprelativeErrorLS;
    relativeErrorGMMOptimal(indexSigma,:,:) = tmprelativeErrorGMMOptimal;
    relativeErrorGMM2Step(indexSigma,:,:) = tmprelativeErrorGMM2Step;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp(['Yo man, im in index....' num2str(indexSigma) '  out of...' num2str(length(sigmaArray))]);
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
end

%% Compare - Relative Error (1 - signal, 2 - rho)
% Ask user if she/he wants to save the figures
toSave = true;

mkdir(savingPath);
save([savingPath,'data.mat']);

ii = 1; % Choose from  1 : 4 (1 - signal Rel. Error, 2 - Rho Rel. Error, 3 - #iter, 4 - CPU Time)
meanrelativeErrorLSIdeal = zeros(length(sigmaArray),1);
meanrelativeErrorOptIdeal = zeros(length(sigmaArray),1);


stdrelativeErrorLSIdeal = zeros(length(sigmaArray),1);
stdrelativeErrorOptIdeal = zeros(length(sigmaArray),1);
proportionsIdeal = zeros(length(sigmaArray),numberRepeats);
%% Compute mean & variance for each sigma
for indexSigma = 1 : length(sigmaArray)
    meanrelativeErrorLSIdeal(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
    meanrelativeErrorOptIdeal(indexSigma) = mean(squeeze(relativeErrorGMMOptimal(indexSigma,:,ii)));

    stdrelativeErrorLSIdeal(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
    stdrelativeErrorOptIdeal(indexSigma) = std(squeeze(relativeErrorGMMOptimal(indexSigma,:,ii)));
    % Ratio
    proportionsIdeal(indexSigma,:) = (squeeze(relativeErrorLS(indexSigma,:,ii))./...
                                squeeze(relativeErrorGMMOptimal(indexSigma,:,ii)));
end

%% 2Steps

ii = 1; % Choose from  1 : 4 (1 - signal Rel. Error, 2 - Rho Rel. Error, 3 - #iter, 4 - CPU Time)
meanrelativeErrorLS2Step = zeros(length(sigmaArray),1);
meanrelativeErrorOpt2Step = zeros(length(sigmaArray),1);


stdrelativeErrorLS2Steps = zeros(length(sigmaArray),1);
stdrelativeErrorOpt2Steps = zeros(length(sigmaArray),1);
proportions2Steps = zeros(length(sigmaArray),numberRepeats);
%% Compute mean & variance for each sigma
for indexSigma = 1 : length(sigmaArray)
    meanrelativeErrorLS2Step(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
    meanrelativeErrorOpt2Step(indexSigma) = mean(squeeze(relativeErrorGMM2Step(indexSigma,:,ii)));

    stdrelativeErrorLS2Steps(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
    stdrelativeErrorOpt2Steps(indexSigma) = std(squeeze(relativeErrorGMM2Step(indexSigma,:,ii)));
    % Ratio
    proportions2Steps(indexSigma,:) = (squeeze(relativeErrorLS(indexSigma,:,ii))./...
                                squeeze(relativeErrorGMM2Step(indexSigma,:,ii)));
end

%% Display Ratio
fig = figure;   
yyaxis left
% [fig] = BoxPlotAsaf(fig, SNR(end:-1:1), proportions(end:-1:1,end:-1:1), 'b*-');
SNRForPlot = SNR(1:end -1);
meanErrror = meanrelativeErrorLS2Step(1:end-1);
meanErrrorGMM = meanrelativeErrorOpt2Step(1:end-1);

proportionsForPlot2Step = proportions2Steps(1:end -1, :);
proportionsForPlotIdeal = proportionsIdeal(1:end -1, :);

[fig,p1,p2] = BoxPlotAsaf(fig, SNRForPlot, proportionsForPlot2Step, 'b*--');
%% Add numerical mean
meansIdeal = mean(proportionsForPlotIdeal,2);
hold on;
yyaxis left

hAx=gca;                                   % retrieve the axes handle
xtk=hAx.XTick;     
yyaxis left

p3 = plot(xtk,meansIdeal, 'ms--');


[indexMeanBigger1] = find(meanrelativeErrorOpt2Step >= 1,1);
hold on;
yyaxis left

plot(1: size(proportionsForPlot2Step,1), ones(size(proportionsForPlot2Step,1),1),'k--');
hold on;
if ~isempty(indexMeanBigger1)
    plot((indexMeanBigger1) * ones(2,1), [max(proportionsForPlot2Step(:)); min(proportionsForPlot2Step(:))],'g');
end

labels = [500,100,10,1,0.1, 0.01];
AA = interp1(SNRForPlot, 1:length(SNRForPlot), labels);
xticks(AA)
xticklabels(labels)
yyaxis left

xlabel('SNR', 'fontsize', 12, 'fontweight','bold');
ylabel('Rel. Error Ratio: LS / GMM' , 'fontsize', 12, 'fontweight','bold');
ylim([0.75,1.75]);
if (ii == 1)
%     title('Ratio Signal Relative Error - LS / GMM');
elseif (ii == 2)
%     title('Ratio \rho Relative Error  - LS / GMM');
elseif (ii == 3)
%     title('Ratio #iter - LS / GMM');
elseif (ii == 4)
%     title('Ratio CPU Time - LS / GMM');
end

yyaxis right
p4 = semilogy(1:length(meanErrror),meanErrror, 'x--');
hold on;
set(gca, 'YScale', 'log')
ylabel('Mean Error of the LS Estimator')

legend([p3,p1, p2, p4],{'Mean - Ideal GMM','Mean - 2-steps GMM', 'Median - 2-steps GMM', 'Mean - rel error LS'}, 'location', 'northwest');


if (toSave)
    if (ii == 1)
        fileName = 'Ratio_Relative_Error_Compare';
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