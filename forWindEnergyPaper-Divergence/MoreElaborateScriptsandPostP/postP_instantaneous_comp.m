close all

for b = 1:length(basePaths)
    figure(Visible="on")
    for c=1:length(cases)
        dataog=windData.(cases{1}).(basePaths{b}){1};
        datahere=windData.(cases{c}).(basePaths{b}){1};
        timearr=0:dt(b):(size(dataog,1)-1)*dt(b);

        plot(timearr,squeeze(datahere(:,1,18,18)),LineWidth=2,Color=colors(c,:))
        hold on
        ylabel('$v (m/s)$',Interpreter='latex')
        box on
        grid on
        grid minor
            grid on
        box on
        
        ax = gca;
        ax.FontSize = 16;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
        xlim([0 (size(dataog,1)-1)*dt(b)])
        legend('$U_{o}$','$U_{div-free}$','$U_{div-free,cons}$','Interpreter','latex')
        xlabel('$t (s)$',Interpreter='latex')
        xlim([0 1000])
   
    end
    
   % keyboard;
    exportgraphics(gcf,['Images\U_',cases{1},'_',basePaths{b},'.png'])

end