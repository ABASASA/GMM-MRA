
load([savingPath 'data.mat']);
toSave = true;
ii = 1; % Choose from  1 : 4 (1 - signal Rel. Error, 2 - Rho Rel. Error, 3 - #iter, 4 - CPU Time)
meanrelativeErrorLS = zeros(length(sigmaArray),1);
meanrelativeErrorOpt = zeros(length(sigmaArray),1);


stdrelativeErrorLS = zeros(length(sigmaArray),1);
stdrelativeErrorOpt = zeros(length(sigmaArray),1);
proportions = zeros(length(sigmaArray),numberRepeats);
%% Compute mean & variance for each sigma
for indexSigma = 1 : length(sigmaArray)
    meanrelativeErrorLS(indexSigma) = mean(squeeze(relativeErrorLS(indexSigma,:,ii)));
    meanrelativeErrorOpt(indexSigma) = mean(squeeze(relativeErrorGMM(indexSigma,:,ii)));

    stdrelativeErrorLS(indexSigma) = std(squeeze(relativeErrorLS(indexSigma,:,ii)));
    stdrelativeErrorOpt(indexSigma) = std(squeeze(relativeErrorGMM(indexSigma,:,ii)));
    % Ratio
    proportions(indexSigma,:) = (squeeze(relativeErrorLS(indexSigma,:,ii))./...
                                squeeze(relativeErrorGMM(indexSigma,:,ii)));
end

fig = figure;
yyaxis left

inxes = 1: length(SNR)-1;
SNRForPlot = SNR(inxes);
meanError = meanrelativeErrorLS(inxes);
proportionsForPlot = proportions(inxes, :);


[fig,p1,p2] = BoxPlotAsaf(fig, SNRForPlot, proportionsForPlot, 'b*--');
% xtick
[indexMeanBigger1] = find(meanrelativeErrorOpt >= 1,1);
hold on;
plot(1: size(proportionsForPlot,1), ones(size(proportionsForPlot,1),1),'k--');
hold on;
if ~isempty(indexMeanBigger1)
    plot((indexMeanBigger1) * ones(2,1), [max(proportionsForPlot(:)); min(proportionsForPlot(:))],'g');
end

labels = [500,100,10,1,0.1, 0.01];
AA = interp1(SNRForPlot, 1:length(SNRForPlot), labels, 'linear','extrap');
xticks(AA)
xticklabels(labels)
xlabel('SNR', 'fontsize', 12, 'fontweight','bold');
ylabel('Ratio of Rel. Error of GMM: 2 / 3 Moments' , 'fontsize', 12, 'fontweight','bold');
ylim([0.9,200]);

yyaxis right
p3 = semilogy(1:length(meanError),meanError, 'x--');
hold on;
set(gca, 'YScale', 'log')
ylabel('Mean Error of the GMM with 2 moments Estimator')

legend([p1,p2,p3],{'Mean', 'Median', 'Mean - rel error LS'}, 'location', 'northwest');


if (ii == 1)
%     title('Ratio Signal Relative Error - LS / GMM');
elseif (ii == 2)
%     title('Ratio \rho Relative Error  - LS / GMM');
elseif (ii == 3)
%     title('Ratio #iter - LS / GMM');
elseif (ii == 4)
%     title('Ratio CPU Time - LS / GMM');
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
