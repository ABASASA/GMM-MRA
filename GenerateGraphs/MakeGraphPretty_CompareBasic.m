filePath = 'E:\Workspace\לימודים\תזה\פורמליזם של הבעיה ל-MRA\GMM_Docs\כתיבת תזה\דמה לסביבת עבודה\MRA-Graphs\2stepGMM_vsLS\LS_With_OptimalGMM_new\';

%% Ideal
load([filePath 'data.mat']);
toSave = true;
savingPath = filePath;
    ii = 1; % Choose from  1 : 4 (1 - signal Rel. Error, 2 - Rho Rel. Error, 3 - #iter, 4 - CPU Time)
    meanrelativeErrorLSIdeal = zeros(length(sigmaArray),1);
    meanrelativeErrorOptIdeal = zeros(length(sigmaArray),1);


    stdrelativeErrorLSIdeal = zeros(length(sigmaArray),1);
    stdrelativeErrorOptIdeal = zeros(length(sigmaArray),1);
    proportionsIdeal = zeros(length(sigmaArray),numberRepeats);
    %% Compute mean & variance for each sigma
    for indexSigma = 1 : length(sigmaArray)
        meanrelativeErrorLSIdeal(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
        meanrelativeErrorOptIdeal(indexSigma) = mean(squeeze(relativeErrorGMM(indexSigma,:,ii)));

        stdrelativeErrorLSIdeal(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
        stdrelativeErrorOptIdeal(indexSigma) = std(squeeze(relativeErrorGMM(indexSigma,:,ii)));
        % Ratio
        proportionsIdeal(indexSigma,:) = (squeeze(relativeErrorLS(indexSigma,:,ii))./...
                                    squeeze(relativeErrorGMM(indexSigma,:,ii)));
    end

    %% 2Steps
    filePath = 'E:\Workspace\לימודים\תזה\פורמליזם של הבעיה ל-MRA\GMM_Docs\כתיבת תזה\דמה לסביבת עבודה\MRA-Graphs\2stepGMM_vsLS\LS_With_2StepGMM-StandardCase\';
    load([filePath 'data.mat']);
    toSave = true;
    savingPath = filePath;

    ii = 1; % Choose from  1 : 4 (1 - signal Rel. Error, 2 - Rho Rel. Error, 3 - #iter, 4 - CPU Time)
    meanrelativeErrorLS2Step = zeros(length(sigmaArray),1);
    meanrelativeErrorOpt2Step = zeros(length(sigmaArray),1);


    stdrelativeErrorLS2Steps = zeros(length(sigmaArray),1);
    stdrelativeErrorOpt2Steps = zeros(length(sigmaArray),1);
    proportions2Steps = zeros(length(sigmaArray),numberRepeats);
    %% Compute mean & variance for each sigma
    for indexSigma = 1 : length(sigmaArray)
        meanrelativeErrorLS2Step(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
        meanrelativeErrorOpt2Step(indexSigma) = mean(squeeze(relativeErrorGMM(indexSigma,:,ii)));

        stdrelativeErrorLS2Steps(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
        stdrelativeErrorOpt2Steps(indexSigma) = std(squeeze(relativeErrorGMM(indexSigma,:,ii)));
        % Ratio
        proportions2Steps(indexSigma,:) = (squeeze(relativeErrorLS(indexSigma,:,ii))./...
                                    squeeze(relativeErrorGMM(indexSigma,:,ii)));
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
% filePath = 'E:\Workspace\לימודים\תזה\פורמליזם של הבעיה ל-MRA\GMM_Docs\כתיבת תזה\דמה לסביבת עבודה\MRA-Graphs\2stepGMM_vsLS\LS_With_OptimalAnalyticalGMM\';
% load([filePath 'data.mat']);
% toSave = false;
% savingPath = filePath;
% ii = 1; % Choose from  1 : 4 (1 - signal Rel. Error, 2 - Rho Rel. Error, 3 - #iter, 4 - CPU Time)
% meanrelativeErrorLS2Step = zeros(length(sigmaArray),1);
% meanrelativeErrorOpt2Step = zeros(length(sigmaArray),1);
% 
% 
% stdrelativeErrorLS2Steps = zeros(length(sigmaArray),1);
% stdrelativeErrorOpt2Steps = zeros(length(sigmaArray),1);
% proportions2Steps = zeros(length(sigmaArray),numberRepeats);
% % %% Compute mean & variance for each sigma
% % for indexSigma = 1 : length(sigmaArray)
% %     meanrelativeErrorLS2Step(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
% %     meanrelativeErrorOpt2Step(indexSigma) = mean(squeeze(relativeErrorGMM(indexSigma,:,ii)));
% % 
% %     stdrelativeErrorLS2Steps(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
% %     stdrelativeErrorOpt2Steps(indexSigma) = std(squeeze(relativeErrorGMM(indexSigma,:,ii)));
% %     % Ratio
% %     proportions2Steps(indexSigma,:) = (squeeze(relativeErrorLS(indexSigma,:,ii))./...
% %                                 squeeze(relativeErrorGMM(indexSigma,:,ii)));
% % end
% 



