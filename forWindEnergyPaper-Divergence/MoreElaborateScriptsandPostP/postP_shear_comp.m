%% SHEAR PROFS
close all
for b = 1:length(basePaths)
    figure(Visible="on")
    for c=1:length(cases)
        subplot(1,2,1)
        plot(meanMean.(cases{c}).(basePaths{b}),[5:10:345]./280,'Color',colors(c,:),LineWidth=2)
        hold on
        box on
        xlabel('$<\overline{u}>$ (m/s)','Interpreter','Latex')
        ylabel('$z/D (-)$',Interpreter='latex')
        ax=gca;
        ax.FontSize=16;
        yticks([0:0.25:1.25])
        ylim([0 1.25])
        yticklabels({'0','0.25','0.5','0.75','1','1.25'})

        box on
        
        ax = gca;
        ax.FontSize = 16;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
        
        yline(170./280,'k--','LineWidth',3)
        yline(170./280-140./280,'k--','LineWidth',3)
        yline(170./280+140./280,'k--','LineWidth',3)
        legend('$U_{o}$','','','','$U_{div-free}$','','','','$U_{div-free,cons}$',Location='northwest',Interpreter='latex')
        subplot(1,2,2)
        plot(meanStd.(cases{c}).(basePaths{b}),[5:10:345]./280,'Color',colors(c,:),LineWidth=2)
        hold on
        xlim([0.1 0.5])
        % grid on
        % grid minor
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

        box on
        
        ax = gca;
        ax.FontSize = 14;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
        
    end
    % keyboard
    exportgraphics(gcf,['Images\MeanstdU_',cases{1},'_',basePaths{b},'.png'])
end