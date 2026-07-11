clear; clc; close all;

% Define seed range and cases
numSeeds = 6;

cases = {'DiffDecay','DiffDecay_divfree'};
% cases = {'SameDecay','SameDecay_divfree'};
% cases = {'SameDecay_2','SameDecay_2_divfree'};
% cases = {'KaimalDecay', 'KaimalDecay_divfree'};
colors=gray(3);
basePaths = {'BR_NABL', 'AR_NABL', 'R_NABL'};  % Multiple base paths

% Initialize structure to store wind data
windData = struct();

% Loop through cases, base paths, and seeds to load data
for c = 1:2%length(cases)
    caseName = cases{c};

    for b = 1:length(basePaths)
        basePath = basePaths{b};

        for i = 1:numSeeds
            seedFolder = ['Seed' num2str(i)];
            filePath = fullfile(caseName, basePath, seedFolder, 'forum_code_wnd_v2.wnd');

            disp(['Checking file: ', filePath]);

            if exist(filePath, 'file') ~= 2
                error(['File not found: ', filePath]);
            end

            windData.(caseName).(basePath){i} = readfile_WND(filePath);
        end
    end
end

% Compute and Plot Mean of Standard Deviations
close all;


% Function to compute mean of standard deviations
% computeMeanStd = @(data) mean(squeeze(std(data(:,1,:,:), 1)), 1);
% computeMeanMean = @(data) mean(squeeze(mean(data(:,1,:,:), 1)), 1);
computeMeanStd = @(data) squeeze(std(data(:,1,:,:), 1));
computeMeanMean = @(data) squeeze(mean(data(:,1,:,:), 1));

computeMeanMeanV = @(data) squeeze(mean(data(:,2,:,:), 1));

get18thRow = @(result) result(18, :);

computeMeanStd = @(data) get18thRow(computeMeanStd(data));
computeMeanMean = @(data) get18thRow(computeMeanMean(data));

computeMeanMeanV = @(data) get18thRow(computeMeanMeanV(data));

% Compute mean std for each case
meanStd = struct();

for c = 1:2%length(cases)
    caseName = cases{c};

    for b = 1:length(basePaths)
        basePath = basePaths{b};

        for i = 1:numSeeds
            meanStd.(caseName).(basePath)(:, i) = computeMeanStd(windData.(caseName).(basePath){i});
            meanMean.(caseName).(basePath)(:, i) = computeMeanMean(windData.(caseName).(basePath){i});
            meanMeanV.(caseName).(basePath)(:, i) = computeMeanMeanV(windData.(caseName).(basePath){i});
        end
        
        % Average across seeds
        meanStd.(caseName).(basePath) = mean(meanStd.(caseName).(basePath), 2);
        meanMean.(caseName).(basePath) = mean(meanMean.(caseName).(basePath), 2);
        meanMeanV.(caseName).(basePath) = mean(meanMeanV.(caseName).(basePath), 2);
    end
end

%%
% SHEAR PROFS
close all
colors=gray(3);
dt=[0.33,0.25,0.25];
Umean=[8.5,14.5,11.5];
for b = 1:length(basePaths)
    figure(Visible="off")
    for c=1:2%length(cases)
        subplot(1,2,1)
        plot(meanMean.(cases{c}).(basePaths{b}),[5:10:345]./280,'Color',colors(c,:),LineWidth=2)
        hold on
        grid on
        grid minor
        box on
        xlabel('$<\overline{u}>$ (m/s)','Interpreter','Latex')
        ylabel('$z/D (-)$',Interpreter='latex')
        ax=gca;
        ax.FontSize=16;
        yticks([0:0.25:1.25])
        ylim([0 1.25])
        yticklabels({'0','0.25','0.5','0.75','1','1.25'})
         grid on
        box on
        
        ax = gca;
        ax.FontSize = 16;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
        
        yline(170./280,'k--','LineWidth',3)
        yline(170./280-140./280,'k--','LineWidth',3)
        yline(170./280+140./280,'k--','LineWidth',3)
        legend('$U_{o}$','','','','$U_{div-free}$',Location='northwest',Interpreter='latex')
        subplot(1,2,2)
        plot(meanStd.(cases{c}).(basePaths{b}),[5:10:345]./280,'Color',colors(c,:),LineWidth=2)
        hold on
        xlim([0.1 0.5])
        grid on
        grid minor
        box on
        xlabel('$\sigma_{U} (m/s)$',Interpreter='latex')
        ylabel('$z/D (-)$',Interpreter='latex')
         grid on
        box on
        
        ax = gca;
        ax.FontSize = 14;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
        
        ax=gca;
        ax.FontSize=16;
        yline(170./280,'k--','LineWidth',3)
        yline(170./280-140./280,'k--','LineWidth',3)
        yline(170./280+140./280,'k--','LineWidth',3)
        yticks([0:0.25:1.25])
        ylim([0 1.25])
        yticklabels({'0','0.25','0.5','0.75','1','1.25'})
         grid on
        box on
        
        ax = gca;
        ax.FontSize = 14;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
        
    end
    % keyboard
    mkdir Images
    exportgraphics(gcf,['Images\MeanstdU_',cases{1},'_',basePaths{b},'.png'])
end