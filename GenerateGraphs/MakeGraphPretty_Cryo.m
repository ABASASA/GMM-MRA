load([savingPath 'data.mat']);
toSave = true;
ii = 1; % Choose from  1 : 4 (1 - signal Rel. Error, 2 - Rho Rel. Error, 3 - #iter, 4 - CPU Time)
meanrelativeErrorLS = zeros(length(ks),1);
meanrelativeErrorOpt = zeros(length(ks),1);


stdrelativeErrorLS = zeros(length(ks),1);
stdrelativeErrorOpt = zeros(length(ks),1);
proportions = zeros(length(ks),numberRepeats);
%% Compute mean & variance for each sigma
for indexK = 1 : length(ks)
    meanrelativeErrorLS(indexK) = mean(squeeze(relativeErrorLS(indexK,:,ii)));
    meanrelativeErrorOpt(indexK) = mean(squeeze(relativeErrorGMM(indexK,:,ii)));

    stdrelativeErrorLS(indexK) = std(squeeze(relativeErrorLS(indexK,:,ii)));
    stdrelativeErrorOpt(indexK) = std(squeeze(relativeErrorGMM(indexK,:,ii)));
    % Ratio
    proportions(indexK,:) = (squeeze(relativeErrorLS(indexK,:,ii))./...
                                squeeze(relativeErrorGMM(indexK,:,ii)));
end


%% Plot histogram - Relative Error

% box plot
 fig = figure;
 yyaxis left

meanError = meanrelativeErrorLS(:);
proportions = proportions(:, :);
[fig,p1,p2] = BoxPlotAsaf(fig, ks, proportions, 'b*-');
[indexMeanBigger1] = find(meanrelativeErrorOpt >= 1,1);
hold on;
plot(1: size(proportions,1), ones(size(proportions,1),1),'k--');
hold on;
if ~isempty(indexMeanBigger1)
    plot((indexMeanBigger1) * ones(2,1), [max(proportions(:)); min(proportions(:))],'g');
end
title('')
xlabel('K', 'fontsize', 12, 'fontweight','bold');
ylabel('Ratio of Rel. Error: LS / GMM' , 'fontsize', 12, 'fontweight','bold');
ylim([0.9,1.5])


yyaxis right
p3 = semilogy(1:length(meanError),meanError, 'x--');
hold on;
set(gca, 'YScale', 'log')
ylabel('Mean Error of the LS Estimator')
ylim([0.0008,0.011])
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